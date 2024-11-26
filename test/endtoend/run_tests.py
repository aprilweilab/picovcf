import unittest
import subprocess
import os

CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(THIS_DIR, "..", "example_vcfs")

def run(cmd):
    print(f"Running: {cmd}")
    return subprocess.check_output(list(map(str, cmd)))

def input(filename):
    return os.path.join(INPUT_DIR, filename)

def quiet_del(filepath):
    if os.path.isfile(filepath):
        os.remove(filepath)


class TestIgdTools(unittest.TestCase):
    def test_unphased_from_igd(self):
        # unphased.example.vcf and msprime.example.vcf are identical, except that the former
        # has been manually "dephased" by changing | to /.
        # So if we convert them both in the correct way, the IGDs should be identical.

        # Convert unphased VCF to IGD
        run(["igdtools", input("unphased.example.vcf"), "--out", "unphased.example.igd"])

        # Convert phased VCF to IGD, then convert phased IGD to unphased IGD
        run(["igdtools", input("msprime.example.vcf"), "--out",  "msprime.example.igd"])
        run(["igdtools", "msprime.example.igd", "--out",  "msprime.example.dephased.igd", "--force-unphased"])

        unphased_result = run(["igdtools", "unphased.example.igd", "--alleles"]).decode("utf-8")
        dephased_result = run(["igdtools", "msprime.example.dephased.igd", "--alleles"]).decode("utf-8")

        self.assertEqual(unphased_result, dephased_result)
        line_count = len(list(filter(lambda x: len(x.strip()) > 0, unphased_result.split("\n"))))
        self.assertEqual(line_count, 8)

        if CLEANUP:
            quiet_del("unphased.example.igd")
            quiet_del("msprime.example.igd")
            quiet_del("msprime.example.dephased.igd")

if __name__ == '__main__':
    unittest.main()
