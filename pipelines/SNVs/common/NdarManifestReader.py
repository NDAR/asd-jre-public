import pandas as pd


class ManifestException(Exception):
    pass

FILE_COLUMN_NAMES = ["data_file%d" % x for x in range(1, 5)]

class NdarSample(object):
    def __init__(self, rows):
        self.rows = rows
        assert len(set(rows["src_subject_id"].values)) == 1
        self.sample_id = rows["src_subject_id"].values[0]
        self.files = []
        for i,r in rows[FILE_COLUMN_NAMES].iteritems():
            self.files.extend(r.dropna().values)
        self.files = list(set(self.files))

    def __repr__(self):
        s = """<NDAR Sample: %s, %d files>""" % (self.sample_id, len(self.files))
        return s

class NdarManifestReader(object):
    def __init__(self, filename):
        df = pd.read_csv(filename, sep="\t", index_col=False)
        df = df[df.index != 0].reindex()
        self.manifest = df
        self.sample_id_col="src_subject_id"
        self.file_counts = {col: self.manifest[col].count() for col in FILE_COLUMN_NAMES}
        self.manifest["family_id"] = map(lambda x: x.split(".")[0], self.manifest[self.sample_id_col])
        self.manifest["rel_id"] = map(lambda x: x.split(".")[1], self.manifest[self.sample_id_col])

    def __repr__(self):
        repr_string = []
        repr_string.append("Samples: %d" % len(self.manifest.src_subject_id.unique()))
        repr_string.append("Total # of files: %d" % sum(self.file_counts.values()))
        for col in FILE_COLUMN_NAMES:
            repr_string.append("  %s: %d" % (col, self.file_counts[col]))
        repr_string.append("Rows: %d" % len(self.manifest))
        return "\n".join(repr_string)

    def get_family(self, family_id):
        mask = self.manifest["family_id"] == family_id
        if sum(mask) == 0:
            raise ManifestException("No samples with family_id %s in manifest!" % family_id)
        sample_list = []
        for sample_id in set(self.manifest[mask][self.sample_id_col].values):
            sample_list.append(self.get_sample(sample_id))
        return sample_list

    def get_sample(self, sample_id):
        mask = self.manifest[self.sample_id_col] == sample_id
        if sum(mask) == 0:
            raise ManifestException("No such sample_id found in %s column of manifest!" % self.sample_id_col)
        else:
            return NdarSample(self.manifest[mask])

    def get_files_by_sample(self, sample_id):
        return self.get_sample(sample_id).files

if __name__ == "__main__":
    m = NdarManifestReader("13835.manifest")
    print m.get_sample("13835.p1")
    print m.get_family("13835")
