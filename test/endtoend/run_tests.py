import unittest
import subprocess
import os

CLEANUP = True

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
INPUT_DIR = os.path.join(THIS_DIR, "..", "example_vcfs")
EXPECT_DIR = os.path.join(THIS_DIR, "expected")

# These tests assume you built igdtools in a build/ directory, or if that fails
# then you can specify `USE_IGDTOOLS_PATH=1` in the environment to force it to use
# the one on the path
IGDTOOLS = (
    "igdtools"
    if (int(os.environ.get("USE_IGDTOOLS_PATH", 0)) == 1)
    else os.path.join(THIS_DIR, "..", "..", "build", "igdtools")
)


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
    def test_merge_igd(self):
        run([IGDTOOLS, input("msprime.example.vcf"), "--out", "msprime.example.igd"])
        run(
            [
                IGDTOOLS,
                "msprime.example.igd",
                "-r",
                "55829-56531",
                "-o",
                "test.merge.pt1.igd",
            ]
        )
        run(
            [
                IGDTOOLS,
                "msprime.example.igd",
                "-r",
                "56554-56812",
                "-o",
                "test.merge.pt2.igd",
            ]
        )
        run(
            [
                IGDTOOLS,
                "test.merge.pt2.igd",
                "--merge",
                "test.merge.pt1.igd",
                "-o",
                "test.merged.igd",
            ]
        )
        result = run([IGDTOOLS, "--alleles", "test.merged.igd"])
        self.assertEqual(
            result.decode("utf-8").split("\n"),
            [
                "POSITION\tREF\tALT\tALT COUNT\tTOTAL",
                "55829\tA\tG\t90\t20000",
                "56531\tA\tG\t329\t20000",
                "56554\tA\tT\t150\t20000",
                "56812\tG\tT\t131\t20000",
                "",
            ],
        )

    def test_subsample_igd(self):
        run([IGDTOOLS, input("msprime.example.vcf"), "--out", "msprime.example.igd"])
        with open("test.subsample.txt", "w") as f:
            f.write("\n".join(["tsk_100", "tsk_6", "tsk_63", "tsk_222"]))
        run(
            [
                IGDTOOLS,
                "-S",
                "test.subsample.txt",
                "msprime.example.igd",
                "-o",
                "test.subsample.igd",
            ]
        )
        result = run([IGDTOOLS, "--individuals", "test.subsample.igd"])
        self.assertEqual(
            result.decode("utf-8").split("\n"),
            ["0: tsk_100", "1: tsk_6", "2: tsk_63", "3: tsk_222", ""],
        )
        result = run([IGDTOOLS, "--alleles", "test.subsample.igd"])
        self.assertEqual(
            result.decode("utf-8").split("\n"),
            [
                "POSITION\tREF\tALT\tALT COUNT\tTOTAL",
                "55829\tA\tG\t0\t8",
                "56531\tA\tG\t1\t8",
                "56554\tA\tT\t1\t8",
                "56812\tG\tT\t0\t8",
                "",
            ],
        )

    def test_unphased_from_igd(self):
        # unphased.example.vcf and msprime.example.vcf are identical, except that the former
        # has been manually "dephased" by changing | to /.
        # So if we convert them both in the correct way, the IGDs should be identical.

        # Convert unphased VCF to IGD
        run([IGDTOOLS, input("unphased.example.vcf"), "--out", "unphased.example.igd"])

        # Convert phased VCF to IGD, then convert phased IGD to unphased IGD
        run([IGDTOOLS, input("msprime.example.vcf"), "--out", "msprime.example.igd"])
        run(
            [
                IGDTOOLS,
                "msprime.example.igd",
                "--out",
                "msprime.example.dephased.igd",
                "--force-unphased",
            ]
        )

        unphased_result = run([IGDTOOLS, "unphased.example.igd", "--alleles"]).decode(
            "utf-8"
        )
        dephased_result = run(
            [IGDTOOLS, "msprime.example.dephased.igd", "--alleles"]
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
                IGDTOOLS,
                input("unphased.example.vcf"),
                "--out",
                "unphased2.example.igd",
            ]
        )

        # Convert phased VCF to IGD, then convert phased IGD to unphased IGD
        run(
            [
                IGDTOOLS,
                input("msprime.example.vcf"),
                "--force-unphased",
                "--out",
                "msprime.example.dephased2.igd",
            ]
        )

        unphased_result = run([IGDTOOLS, "unphased2.example.igd", "--alleles"]).decode(
            "utf-8"
        )
        dephased_result = run(
            [IGDTOOLS, "msprime.example.dephased2.igd", "--alleles"]
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
        run([IGDTOOLS, input("multi.vcf"), "-o", "multi.igd"])
        run(
            [
                IGDTOOLS,
                "multi.igd",
                "--drop-multi-sites",
                "-o",
                "multi.dropped_sites.igd",
            ]
        )

        igd_tools_result = run(
            [IGDTOOLS, "multi.dropped_sites.igd", "--alleles"]
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
        run([IGDTOOLS, input("multi.vcf"), "-o", "multi.igd"])
        run(
            [
                IGDTOOLS,
                "multi.igd",
                "--drop-non-snv-sites",
                "-o",
                "multi.only_snv_sites.igd",
            ]
        )

        igd_tools_result = run(
            [IGDTOOLS, "multi.only_snv_sites.igd", "--alleles"]
        ).decode("utf-8")
        with open(expect("multi.only_snv_sites.af.txt"), "rb") as fexpect:
            expected_result = fexpect.read().decode("utf-8")
        self.assertEqual(igd_tools_result, expected_result)

        # VARIANTS
        run([IGDTOOLS, "multi.igd", "--drop-non-snvs", "-o", "multi.only_snvs.igd"])

        igd_tools_result = run([IGDTOOLS, "multi.only_snvs.igd", "--alleles"]).decode(
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
