#!/usr/bin/env python3
import hashlib
from os import listdir
from os.path import join, dirname

from pytest import fixture

from bin.count_reads_per_region import QC, Counters


@fixture
def design_file(tmp_path):
    """Making this a fixture so that the test is consistent if the original file ever changes"""

    design_file = tmp_path / "PFA_GRC1_v1.0.annotation.regions.txt"
    design_file.write_text(
        "Pf3D7_01_v3_145515_294I_A:1-205\n"
        "Pf3D7_01_v3_180554_D714N:1-245\n"
        "Pf3D7_01_v3_535211_2521F_A:1-231\n"
        "Pf3D7_02_v3_470013_G75E_A:1-246\n"
        "Pf3D7_02_v3_714480_D258G:1-290\n"
        "Pf3D7_03_v3_155697_150P_B:1-221\n"
        "Pf3D7_03_v3_656861_129V_A:1-234\n"
        "Pf3D7_04_v3_139051_K438N_A:1-193\n"
        "Pf3D7_04_v3_426436_D560A:1-214\n"
        "Pf3D7_04_v3_531138_A992E_A:1-239\n"
        "DHFR_16_51_59:1-245\n"
        "DHFR_306:1-193\n"
        "Pf3D7_04_v3_881571_1081R_A:1-200\n"
        "Pf3D7_04_v3_1037656_2776I:1-225\n"
        "Pf3D7_05_v3_172801_E218K_A:1-248\n"
        "Pf3D7_05_v3_369740_907L_A:1-200\n"
        "MDR1_86:1-248\n"
        "Pf3D7_05_v3_1204155_1338I_A:1-228\n"
        "Pf3D7_06_v3_900278_P696S_A:1-224\n"
        "Pf3D7_06_v3_1289212_125T_A:1-196\n"
        "CRT_72_74_75_76:1-205\n"
        "CRT_220:1-249\n"
        "CRT_326:1-244\n"
        "CRT_371:1-213\n"
        "Pf3D7_07_v3_619957_675R_A:1-201\n"
        "Pf3D7_07_v3_704373_389E:1-191\n"
        "Pf3D7_07_v3_1066698_G483S_A:1-205\n"
        "Pf3D7_07_v3_1256331_L321F_A:1-202\n"
        "Pf3D7_07_v3_1358910:1-240\n"
        "Pf3D7_08_v3_150033_1315I_A:1-244\n"
        "Pf3D7_08_v3_399774_421K_A:1-195\n"
        "Pf3D7_08_v3_417335_R244K_A:1-198\n"
        "DHPS_436_437:1-218\n"
        "Pf3D7_08_v3_549993:1-203\n"
        "Pf3D7_08_v3_1056829_L474I_A:1-232\n"
        "Pf3D7_08_v3_1314831_1342K_A:1-226\n"
        "Pf3D7_09_v3_452690_1018I:1-201\n"
        "Pf3D7_09_v3_900277_1534E_A:1-199\n"
        "Pf3D7_10_v3_361684:1-190\n"
        "Pf3D7_10_v3_1383789_N114H:1-241\n"
        "Pf3D7_10_v3_1386850_927K_A:1-190\n"
        "Pf3D7_11_v3_477922_H147Y_A:1-191\n"
        "Pf3D7_11_v3_1006911_D124E_B:1-241\n"
        "Pf3D7_11_v3_1020397_G700E_A:1-228\n"
        "Pf3D7_11_v3_1295068_E405K:1-221\n"
        "Pf3D7_11_v3_1815412_E765Q:1-248\n"
        "Pf3D7_11_v3_1935031_I139L:1-197\n"
        "Pf3D7_12_v3_858501_Q469K:1-195\n"
        "Pf3D7_12_v3_974663:1-202\n"
        "Pf3D7_12_v3_1667593_2381N_A:1-210\n"
        "Pf3D7_12_v3_2171901_V140D_A:1-208\n"
        "Pf3D7_13_v3_388365_S1236R:1-195\n"
        "Pf3D7_13_v3_1056452_1234D:1-221\n"
        "Pf3D7_13_v3_1419519:1-241\n"
        "K13_resistance_1:1-245\n"
        "K13_resistance_3:1-215\n"
        "K13_resistance_5:1-249\n"
        "Pf3D7_13_v3_1867630_M4911I:1-229\n"
        "Pf3D7_13_v3_2377887_2002S_A:1-201\n"
        "Pf3D7_13_v3_2573828_I1153M:1-225\n"
        "Pf3D7_14_v3_137622_1179V_A:1-246\n"
        "PlasIV_ref:1-260\n"
        "PlasII_ref:1-260\n"
        "Pf3D7_14_v3_1757603_D1365G:1-248\n"
        "Pf3D7_14_v3_2164225_2830S_B:1-234\n"
        "Pf3D7_14_v3_2733656_557C_A:1-240\n"
        "Pf3D7_14_v3_3126219:1-232\n"
        "Pf_PF3D7_1460900-1_Pf3D7_14:1-213"
    )
    yield design_file


@fixture
def plex_file(tmp_path):
    plex_file = tmp_path / "201119_GRC1.plex"
    plex_file.write_text("201119,1,64\n" "201119,1,71")
    yield plex_file


@fixture
def tmp_out(tmp_path):
    count_reads_per_region = tmp_path / "count_reads_per_region.txt"
    yield count_reads_per_region


MAPQ = 10
CRAM_PATHS = join(dirname(__file__), "crams")
CRAM_FILE_1 = join(CRAM_PATHS, "201119_1#64.cram")
CRAM_FILE_2 = join(CRAM_PATHS, "201119_1#71.cram")


def test_counters():
    c = Counters(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
    assert c.zero == 0.1
    assert c.ten == 0.2
    assert c.ten_f == 0.3
    assert c.ten_r == 0.4
    assert c.ten_frag == 0.5
    assert c.ten_fragboth == 0.6


def test_qc_read_text_file(design_file):
    lines = QC._read_text_file(design_file)
    assert len(lines) == 68


def test_qc_read_csv_file(plex_file):
    lines = QC._read_csv_file(plex_file)
    assert lines[0][0] == "201119"
    assert lines[1][2] == "71"
    assert len(lines) == 2
    assert len(lines[0]) == 3
    assert len(lines[1]) == 3


def test_qc_run_samtools(design_file, plex_file, tmp_out):
    with QC(design_file, plex_file, CRAM_PATHS, MAPQ, tmp_out) as qc:
        lines = qc._run_samtools(f"-c -F0xB04 {CRAM_FILE_1}")
    assert lines[0] == "1178"


def test_qc_start_task(design_file, plex_file, tmp_out):
    with QC(design_file, plex_file, CRAM_PATHS, MAPQ, tmp_out) as qc:
        task = qc._start_task(["-c", "-F0xB04", CRAM_FILE_1])
        lines = task.get()
    assert lines[0] == "1178"


def test_calculate_samtools_count(design_file, plex_file, tmp_out):
    with QC(design_file, plex_file, CRAM_PATHS, MAPQ, tmp_out) as qc:
        lines = qc._run_samtools(CRAM_FILE_1 + " K13_resistance_3")
        result = qc._calculate_samtools_count(lines, "194", "180")
        assert result == 6


def test_get_samtools_count(design_file, plex_file, tmp_out):
    with QC(design_file, plex_file, CRAM_PATHS, MAPQ, tmp_out) as qc:
        assert qc._get_samtools_count([]) == 0
        assert qc._get_samtools_count(["0"]) == 0
        assert qc._get_samtools_count(["12345", "67890"]) == 12345


def test_calculate_samtools_frag_counts(design_file, plex_file, tmp_out):
    with QC(design_file, plex_file, CRAM_PATHS, MAPQ, tmp_out) as qc:
        lines = qc._run_samtools(CRAM_FILE_1 + " K13_resistance_3")
        result = qc._calculate_samtools_frag_counts(lines)
        assert result == (4, 4)


def test_run(design_file, plex_file, tmp_out):
    with QC(design_file, plex_file, CRAM_PATHS, MAPQ, tmp_out) as qc:
        qc.run()
    with open(tmp_out, "rb") as fh:
        md5 = hashlib.md5(fh.read()).hexdigest()
    assert md5 == "50dd0ee2f8000458a5b1211f3fdbba89"

