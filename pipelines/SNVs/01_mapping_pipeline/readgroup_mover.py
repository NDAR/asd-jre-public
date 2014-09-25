import pysam
import argparse
import cPickle
import sys
import time

parser = argparse.ArgumentParser("A program to re-attach *multiple* readgroups from a bam file to a sam file")
parser.add_argument("command")
parser.add_argument("--prefix", help="Prefix all RG tags (e.g., with filename or sampleID); can help make non-unique RGs unique.", default=None)
parser.add_argument("--dictfile")
parser.add_argument("--lookups")
parser.add_argument("--input", default='stdin')
parser.add_argument("--output", default='stdout')
#parser.add_argument("--use-trie", default=False, )

args = parser.parse_args()

HEADER_FIELDS = ["ID", "PL", "PU", "LB", "DS", "DT", "SM", "CN"]

if args.command == 'create':
    bam_file = pysam.Samfile(args.input,'rb',check_header=False, check_sq=False)
    rg_header = bam_file.header["RG"]
    rg_dict = {}
    rg_lookup = {}
    rg_key = 0
    t_start = time.time()
    t1 = t_start
    for ix, read in enumerate(bam_file.fetch(until_eof=True)):
        if ix % 1e6 == 0:
            t_now = time.time()
            print "[INFO]\tReadgroup_mover.py\tProcessed {read_count} reads\tIteration time: {iter_time}\tTotal time: {total_time}"\
                .format(read_count=ix,
                        iter_time=t_now-t1,
                        total_time=t_now-t_start)
            t1 = t_now
        try:
            rg = rg_lookup[dict(read.tags)["RG"]]
        except KeyError:
            rg_lookup[dict(read.tags)["RG"]] = rg_key
            rg = rg_key
            rg_key += 1
        rg_dict[read.qname] = rg
    t1 = time.time()

    if args.prefix:
        # prefix the ID in the header with the --prefix line
        for i in rg_header:
            i["ID"] = "%s.%s" % (args.prefix,i["ID"])
        # fix the rg_lookup
        new_lookup_keys = ["%s.%s" % (args.prefix, i) for i in rg_lookup.keys()]
        rg_lookup = dict(zip(new_lookup_keys, rg_lookup.values()))


    cPickle.dump(rg_dict,open(args.output,'w'), 2)
    cPickle.dump(rg_lookup,open(args.output + ".lookups",'w'), 2)
    cPickle.dump(rg_header,open(args.output + ".header",'w'), 2)
    print "Succesfully wrote %d reads to %s. Write time: %d" % (ix, args.output, t1-time.time())
    print "Done with readgroup create stage. Total time: %d" % (t_start - time.time())

elif args.command == 'translate':
    rg_dict = cPickle.load(open(args.dictfile,'r'))
    rg_lookup = cPickle.load(open(args.dictfile+".lookups",'r'))
    rg_lookup = dict(zip(rg_lookup.values(), rg_lookup.keys()))
    rg_header = cPickle.load(open(args.dictfile+".header",'r'))
    if args.input == 'stdin':
        line = sys.stdin.readline()
        while line:
            # check if a header
            if line[0] == "@":
                if line.startswith("@RG\t"):
                    #skip existing RG lines, these will be created
                    continue
                else:
                    sys.stdout.write(line)
            else:
                # write new RG lines:
                for h in rg_header:
                    # order common keys for @RG line
                    out_line = "@RG\t"
                    for field in HEADER_FIELDS:
                        if field in h:
                            out_line += "%s:%s\t" % (field, h[field])
                            del h[field] # remove this from dict so as to mark it as added to line
                    # get the remaining keys
                    for field in h.keys():
                        out_line += "%s:%s\t" % (field, h[field])
                    # write out the @RG line, removing trailing \t
                    sys.stdout.write(out_line.rstrip("\t") + "\n")
                break
            line = sys.stdin.readline()

        # now keep runnign through rest of file
        while line:
            rg = rg_lookup[rg_dict[line.split("\t",1)[0]]]
            sys.stdout.write("%s\tRG:Z:%s\n" % (line.rstrip("\n"), rg))
            line = sys.stdin.readline()
    else:
        print "Not implemented"

else:
    print "Error, command not understand. Please use either create or translate commands"
    sys.exit(1)
