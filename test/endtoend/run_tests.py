import unittest
import subprocess
import os

CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(THIS_DIR, "..", "example_vcfs")
EXPECT_DIR = os.path.join(THIS_DIR, "expected")


def run(cmd):
    print(f"Running: {cmd}")
    return subprocess.check_output(list(map(str, cmd)))


def input(filename):
    return os.path.join(INPUT_DIR, filename)


def expect(filename):
    return os.path.join(EXPECT_DIR, filename)


def quiet_del(filepath):
    if os.path.isfile(filepath):
        os.remove(filepath)


class TestIgdTools(unittest.TestCase):
    def test_unphased_from_igd(self):
        # unphased.example.vcf and msprime.example.vcf are identical, except that the former
        # has been manually "dephased" by changing | to /.
        # So if we convert them both in the correct way, the IGDs should be identical.

        # Convert unphased VCF to IGD
        run(
            ["igdtools", input("unphased.example.vcf"), "--out", "unphased.example.igd"]
        )

        # Convert phased VCF to IGD, then convert phased IGD to unphased IGD
        run(["igdtools", input("msprime.example.vcf"), "--out", "msprime.example.igd"])
        run(
            [
                "igdtools",
                "msprime.example.igd",
                "--out",
                "msprime.example.dephased.igd",
                "--force-unphased",
            ]
        )

        unphased_result = run(["igdtools", "unphased.example.igd", "--alleles"]).decode(
            "utf-8"
        )
        dephased_result = run(
            ["igdtools", "msprime.example.dephased.igd", "--alleles"]
        ).decode("utf-8")

        self.assertEqual(unphased_result, dephased_result)
        line_count = len(
            list(filter(lambda x: len(x.strip()) > 0, unphased_result.split("\n")))
        )
        self.assertEqual(line_count, 8)

        if CLEANUP:
            quiet_del("unphased.example.igd")
            quiet_del("msprime.example.igd")
            quiet_del("msprime.example.dephased.igd")

    def test_unphased_from_vcf(self):
        # Same as test_unphased_from_igd, except we unphase the file while converting from VCF

        # Convert unphased VCF to IGD
        run(
            [
                "igdtools",
                input("unphased.example.vcf"),
                "--out",
                "unphased2.example.igd",
            ]
        )

        # Convert phased VCF to IGD, then convert phased IGD to unphased IGD
        run(
            [
                "igdtools",
                input("msprime.example.vcf"),
                "--force-unphased",
                "--out",
                "msprime.example.dephased2.igd",
            ]
        )

        unphased_result = run(
            ["igdtools", "unphased2.example.igd", "--alleles"]
        ).decode("utf-8")
        dephased_result = run(
            ["igdtools", "msprime.example.dephased2.igd", "--alleles"]
        ).decode("utf-8")

        self.assertEqual(unphased_result, dephased_result)
        line_count = len(
            list(filter(lambda x: len(x.strip()) > 0, unphased_result.split("\n")))
        )
        self.assertEqual(line_count, 8)

        if CLEANUP:
            quiet_del("unphased2.example.igd")
            quiet_del("msprime.example.dephased2.igd")

    def test_multi_allelic(self):
        # Verify that removing multi-allelic sites works
        run(["igdtools", input("multi.vcf"), "-o", "multi.igd"])
        run(
            [
                "igdtools",
                "multi.igd",
                "--drop-multi-sites",
                "-o",
                "multi.dropped_sites.igd",
            ]
        )

        igd_tools_result = run(
            ["igdtools", "multi.dropped_sites.igd", "--alleles"]
        ).decode("utf-8")
        with open(expect("multi.dropped_sites.af.txt"), "rb") as fexpect:
            expected_result = fexpect.read().decode("utf-8")

        self.assertEqual(igd_tools_result, expected_result)

        if CLEANUP:
            quiet_del("multi.igd")
            quiet_del("multi.dropped_sites.igd")

    def test_snv_filters(self):
        # Verify that removing non-snvs works for both sites and variants.

        # SITES
        run(["igdtools", input("multi.vcf"), "-o", "multi.igd"])
        run(
            [
                "igdtools",
                "multi.igd",
                "--drop-non-snv-sites",
                "-o",
                "multi.only_snv_sites.igd",
            ]
        )

        igd_tools_result = run(
            ["igdtools", "multi.only_snv_sites.igd", "--alleles"]
        ).decode("utf-8")
        with open(expect("multi.only_snv_sites.af.txt"), "rb") as fexpect:
            expected_result = fexpect.read().decode("utf-8")
        self.assertEqual(igd_tools_result, expected_result)

        # VARIANTS
        run(["igdtools", "multi.igd", "--drop-non-snvs", "-o", "multi.only_snvs.igd"])

        igd_tools_result = run(["igdtools", "multi.only_snvs.igd", "--alleles"]).decode(
            "utf-8"
        )
        with open(expect("multi.only_snvs.af.txt"), "rb") as fexpect:
            expected_result = fexpect.read().decode("utf-8")
        self.assertEqual(igd_tools_result, expected_result)

        if CLEANUP:
            quiet_del("multi.igd")
            quiet_del("multi.only_snv_sites.igd")
            quiet_del("multi.only_snvs.igd")


if __name__ == "__main__":
    unittest.main()
