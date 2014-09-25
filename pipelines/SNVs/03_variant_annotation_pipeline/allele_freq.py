import redis
import pandas as pd
from cStringIO import StringIO
import argparse 
import sys
import snphwe

pd.set_option("max_columns", -1)
pd.set_option("line_width", 5000)
pd.set_option("max_colwidth", 5000)
pd.set_option("max_rows", 2000)

VERSION = "0.1"

def read_vcf(vcf_filename):
    s = StringIO()
    vcf_header_lines = ""
    with open(vcf_filename) as f:
        for line in f:
            if line[0:2] == "##":
                vcf_header_lines += line
            elif line[0] == '#':
                columns = line.split()
            else:
                s.write("%s\n" % line.strip("\t\n"))
    s.seek(0)
    #df = pd.read_csv(s, sep="\t",names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","DATA"])
    df = pd.read_csv(s, sep="\t",names=columns)
    return df, vcf_header_lines, columns

def write_vcf(filename, df, header_lines, columns):
    with open(filename, 'w') as out_f:
        out_f.writelines(header_lines)
        df[columns].to_csv(out_f, sep="\t", index=False)

def missing_gt(row, samples):
    for s in samples:
        if row[s][0:2] == "NA":
            return True
    return False


def parse_samples(samples):
    sample_dict = {}
    for s in samples:
        fID, rel = s.split(".")
        if rel == "fa":
            sample_dict["fa"] = s
        elif rel == "mo":
            sample_dict["mo"] = s
        elif rel[0] == "p":
            sample_dict["pro"] = s
        elif rel[0] == "s":
            sample_dict["sib"] = s
        else:
            Exception("Unknown sample ID encountered!")
    return sample_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Frequency_stats.py v%s" % VERSION)
    parser.add_argument("--redis-host", default="localhost")
    parser.add_argument("--namespace", default="")
    parser.add_argument("--verbose", default=False, action="store_true")
    subparsers = parser.add_subparsers(help='commands', dest="command")

    count_parser = subparsers.add_parser('count', help='Ingest a VCF into the redis DB')
    count_parser.add_argument("vcf_file")

    #hwe_parser = subparsers.add_parser('hwe', help='Calculate HWE p-values for the redis DB')
    #hwe_parser.add_argument("--min-freq", default=0.01)

    annotate_parser = subparsers.add_parser('annotate', help='Annotate a VCF with AC info.')
    annotate_parser.add_argument("--hwe-min-freq", default=0.01)
    annotate_parser.add_argument("--force-hwe-calculation", default=False, action="store_true")
    annotate_parser.add_argument("vcf_file")
    annotate_parser.add_argument("vcf_out")
    annotate_parser.add_argument("n_alleles", type=int)
    
    args = parser.parse_args()

    # read VCF
    vcf, header, columns = read_vcf(args.vcf_file)
    samples = parse_samples(columns[9:])
    
    # filter VCF:
    #vcf["is_mendelian"] = map(lambda x: "Mendelian=True" in x[1]["INFO"], vcf.iterrows())
    #mask = vcf.is_mendelian == True

    # set up redis
    r_server = redis.Redis(args.redis_host)

    if args.command == "count":
        for ix, row in vcf.iterrows():
            hash_key = "%s:%s:%s:%s/%s" % (args.namespace, row["#CHROM"], row["POS"], row["REF"], row["ALT"])

            fmt = row["FORMAT"].split(":")
            fa = dict(zip(fmt, row[samples["fa"]].split(":")))
            mo = dict(zip(fmt, row[samples["mo"]].split(":")))
            
            if "GT" in fmt:
                if fa["GT"] != "NA":
                    if ("FGT" in fa) and (fa["FGT"] != "NA"):
                        r_server.hincrby(hash_key, fa["FGT"])
                    else:
                        r_server.hincrby(hash_key, fa["GT"])
                if mo["GT"] != "NA":
                    if ("FGT" in mo) and (mo["FGT"] != "NA"):
                        r_server.hincrby(hash_key, mo["FGT"])
                    else:
                        r_server.hincrby(hash_key, mo["GT"])

    elif args.command == "annotate":
        out_str = []
        for ix, row in vcf.iterrows():
            hash_key = "%s:%s:%s:%s/%s" % (args.namespace, row["#CHROM"], row["POS"], row["REF"], row["ALT"])
            data = r_server.hgetall(hash_key)
            ref_total_ac = 0
            alt_total_ac = 0
            hwe_p = None
            obs_hets = 0
            obs_hom2 = 0
            obs_other = 0
            for key, value in data.iteritems():
                if key == "hwe":
                    if not args.force_hwe_calculation:
                        hwe_p = value
                else:
                    cnt = int(value)
                    if (key == "0/1"):
                        ref_total_ac += cnt
                        alt_total_ac += cnt
                        obs_hets = cnt
                    elif (key == "1/1"):
                        alt_total_ac += (2 * cnt)
                        obs_hom2 = cnt
                    elif (key == "0/0"):
                        ref_total_ac += (2 * cnt)
                    else:
                        obs_other = cnt
                        if args.verbose:
                            sys.stderr.write("unexpected gt at %s\n" % hash_key)

            obs_hom1 = int((args.n_alleles/2.) - obs_hets - obs_hom2 - obs_other)

            alt_af = float(alt_total_ac)/args.n_alleles
            ref_af = 1 - alt_af
            maf = min(alt_af, ref_af)

            if (hwe_p == None) and (maf > args.hwe_min_freq ):
                hwe_p = snphwe.SNPHWE(obs_hets, obs_hom1, obs_hom2)            
            else:
                hwe_p = "NC"
            
            _ = r_server.hset(hash_key, "hwe", hwe_p)                
            #if alt_total_ac > 0:
            
            out_str.append(";POP_GT_COUNTS=%d,%d,%d;POP_ALT_AC=%d;POP_HWE_p=%s" % (obs_hom1, obs_hets, obs_hom2, alt_total_ac, str(hwe_p)))
            #else:
            #    out_str.append("")

        vcf["INFO"] += out_str
        
        header += '##allele_freq.py=<Version=%s Nik Krumm 2014 nkrumm@uw.edu,Command="%s">\n' % (VERSION, " ".join(sys.argv))
        header += '##allele_freq.py=<n_alleles=%d,population=%d>\n' % (args.n_alleles,args.n_alleles/2)
        header += '##INFO=<ID=POP_GT_COUNTS,Number=3,Type=Integer,Description="Population counts for 0/0, 0/1 and 1/1 genotypes; 0/0 inferred from n_alleles parameter">\n'
        header += '##INFO=<ID=POP_ALT_AC,Number=1,Type=Integer,Description="Population alternate allele count">\n'
        header += '##INFO=<ID=POP_HWE_p,Number=1,Type=Float,Description="Hardy-Weinberg Equilibrium probability (code based on PMID:15789306)">\n'
        write_vcf(args.vcf_out, vcf, header, columns)