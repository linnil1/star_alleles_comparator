import subprocess

import pytest
from star_allele_comp import *
from star_allele_comp.command import entrypoint


def test_allele():
    # any string is OK
    HlaAllele("A*")
    assert HlaAllele("A*").resolution == 0
    HlaAllele("A*new")
    assert HlaAllele("A*new").resolution == 1
    HlaAllele("A*01:02:03:04")
    assert HlaAllele("A*01:02:03:04").resolution == 4
    # Only gene name is not OK
    with pytest.raises(AlleleError):
        HlaAllele("A")
    # So does : in HLA
    with pytest.raises(AlleleError):
        HlaAllele("A*123::gg")
    with pytest.raises(AlleleError):
        HlaAllele("A*123:01:")
    HlaAllele("A*123#non-digit-char-is-treated-as-comment.1.2.3")

    # test ignore_suffix
    HlaAllele("A*123:01w")
    assert HlaAllele("A*123:01w", ignore_suffix=True).trim_resolution(2) == HlaAllele("A*123:01")
    assert HlaAllele("A*123:01", ignore_suffix=True).trim_resolution(2) == HlaAllele("A*123:01")
    KirAllele("2DL1*1230101w")
    assert KirAllele("2DL1*1230101w", ignore_suffix=True).trim_resolution(3) == KirAllele("2DL1*1230101")
    assert KirAllele("2DL1*1230101", ignore_suffix=True).trim_resolution(3) == KirAllele("2DL1*1230101")

    # KIR
    KirAllele("KIR2DL1*")
    assert KirAllele("KIR2DL1*").resolution == 0
    KirAllele("KIR2DL1*new")
    assert KirAllele("KIR2DL1*new").resolution == 1
    KirAllele("KIR2DL1*001")
    KirAllele("KIR2DL1*00105")
    KirAllele("KIR2DL1*00105N")
    KirAllele("KIR2DL1*00105NN")
    KirAllele("KIR2DL1*0010101")
    assert KirAllele("KIR2DL1*0010203").resolution == 3
    KirAllele("KIR2DL1*00101001#ok-comment")
    with pytest.raises(AlleleError):
        KirAllele("KIR2DL1*0010")
    with pytest.raises(AlleleError):
        KirAllele("KIR2DL1*0010N")

    assert KirAllele("KIR2DL1*0010203").trim_resolution(2) == KirAllele("KIR2DL1*00102")
    assert KirAllele("KIR2DL1*00102").trim_resolution(2) == KirAllele("KIR2DL1*00102")
    assert KirAllele("KIR2DL1*001").trim_resolution(2) == KirAllele("KIR2DL1*001")


def test_compare_cohort():
    cohort1 = {
        "id": [
            "A*01:02:03:04",
            "B*01:02:03:05",
            "C*01:02:03",
            "D*01:02:03:04",
            "E*03:02:03:04",
            "G*01:02:03:04",
            "H*01:02:03:04",
        ]
    }
    cohort2 = {
        "id": [
            "A*01:02:03:04",
            "B*01:02:03:04",
            "C*01:02:03:04",
            "D*01:02:03",
            "E*01:02:03:04",
            "F*01:02:03:04",
            "H*01:02:06:04",
        ]
    }
    result = compare_cohort(cohort1, cohort2, "hla")
    assert [i.match_str for i in result["id"]] == [
        "4",
        "3",
        "3",
        "3",
        "0",
        "FP",
        "FN",
        "2",
    ]

    cohort2 = {"id1": cohort2["id"]}
    result = compare_cohort(cohort1, cohort2, "hla")
    assert all(i.match_str == "FN" for i in result["id"])
    assert all(i.match_str == "FP" for i in result["id1"])

    cohort1 = {
        "id": [
            "A*01:02:03:04",
            "A*01:02:03:05",
            "A*08:02:03:06",
            "B*01:02:03:04",
            "B*02:02:03:04",
        ]
    }
    cohort2 = {
        "id": ["A*01:02:02:04", "A*01:02:03:04", "A*08:02:03:04", "B*02:02:03:08"]
    }
    result = compare_cohort(cohort1, cohort2, "HLA")
    assert set([i.match_str for i in result["id"] if i.gene == "A"]) == set(
        ["4", "2", "3"]
    )
    assert set([i.match_str for i in result["id"] if i.gene == "B"]) == set(["FN", "3"])


def test_compare_method():
    method_cohort = {
        "method1": {
            "id": [
                "2DL1*0010203",
                "2DL2*00102",
                "2DL3*0010203",
                "2DL4*0010203",
                "2DL5*0010203",
            ]
        },
        "method2": {
            "id": [
                "2DL1*0010208",
                "2DL2*00102",
                "2DL3*0020203",
                "2DL4*001",
                "2DL5*0010203",
            ]
        },
        "method3": {"id": ["2DL1*0010208", "2DL2*00108", "2DL3*0030203", "2DL5A*001"]},
    }
    result = compare_method(method_cohort, "method1", "kir")
    assert [i.match_res for i in result["method1"]["id"]] == [3, 2, 3, 3, 3]
    assert [i.match_str for i in result["method2"]["id"]] == ["2", "2", "0", "1", "3"]
    assert [i.match_str for i in result["method3"]["id"]] == [
        "2",
        "1",
        "0",
        "FN",
        "FN",
        "FP",
    ]


def test_summary():
    method_cohort = {
        "method1": {
            "id": [
                "2DL1*0010203",
                "2DL2*00102",
                "2DL3*0010203",
                "2DL4*0010203",
                "2DL5*0010203",
            ]
        },
        "method2": {
            "id": [
                "2DL1*0010208",
                "2DL2*00102",
                "2DL3*0020203",
                "2DL4*001",
                "2DL5*0010203",
            ]
        },
        "method3": {"id": ["2DL1*0010208", "2DL2*00108", "2DL3*0030203", "2DL5A*001"]},
    }
    result = compare_method(method_cohort, "method1", "kir")
    result = result.to_dataframe()
    result_df = table_summarize(result)
    result_df = result_df.reset_index()
    print(result_df)

    def get_value(metric, method_want, resolution_want):
        return result_df.query(
            f"method == '{method_want}' and Resolution == '{resolution_want}'"
        ).iloc[0][metric]

    # [i.match_res for i in result["method1"]["id"]] == [3, 2, 3, 3, 3]
    assert get_value("Accuracy", "method1", 3) == 1.0
    assert get_value("Accuracy", "method1", "FN") == 0.0

    # [i.match_str for i in result["method2"]["id"]] == ["2", "2", "0", "1", "3"]
    assert get_value("num_match", "method2", "3") == 1  # 3
    assert get_value("num_match", "method2", "2") == 3  # 3 2 2
    assert get_value("num_match", "method2", "1") == 4  # 3 2 2 1
    assert get_value("num_match", "method2", "0") == 5  # 3 2 2 1 0
    assert get_value("Accuracy", "method2", "3") == 1 / 4
    assert get_value("Accuracy", "method2", "2") == 3 / 5
    assert get_value("Accuracy", "method2", "1") == 4 / 5

    # [i.match_str for i in result["method3"]["id"]] == ["2", "1", "0", "FN", "FN", "FP"]
    assert get_value("num_match", "method3", "FN") == 2
    assert get_value("num_match", "method3", "FP") == 1
    assert get_value("num_match", "method3", "0") == 3  # 2 1 0


@pytest.fixture
def cohort():
    cohort1 = {
        "id": [
            "A*01:02:03:04",
            "B*01:02:03:05",
            "C*01:02:03",
            "D*01:02:03:04",
            "E*03:02:03:04",
            "G*01:02:03:04",
            "H*01:02:03:04",
        ]
    }
    cohort2 = {
        "id": [
            "A*01:02:03:04",
            "B*01:02:03:04",
            "C*01:02:03:04",
            "D*01:02:03",
            "E*01:02:03:04",
            "F*01:02:03:04",
            "H*01:02:06:04",
        ]
    }
    return {
        "method1": cohort1,
        "method2": cohort2,
    }


def test_summary_no_error(cohort, tmp_path):
    result1 = compare_cohort(cohort["method1"], cohort["method2"], "hla")
    result1_df = result1.to_dataframe()
    print(result1)
    print_all_summary(result1_df, include_raw=True)
    save_all_summary(result1_df, tmp_path)

    # cohort
    result2 = compare_method(cohort, "method1", "hla")
    result2_df = result2.to_dataframe()
    print_all_summary(result2_df, include_raw=True)
    save_all_summary(result2_df, tmp_path)


def test_plot_no_error(cohort, tmp_path):
    result1 = compare_method(cohort, "method1", "hla")
    result1_df = result1.to_dataframe()
    for fig in list(plot_summary(result1_df)):
        fig.write_image(tmp_path / "fig1.png")

    result2 = compare_cohort(cohort["method1"], cohort["method2"], "hla")
    result2_df = result2.to_dataframe()
    for fig in list(plot_summary(result2_df)):
        fig.write_image(tmp_path / "fig2.png")


def test_command_no_error(tmp_path):
    command = f"star_allele_comp tests/test1.csv tests/test2.csv --family hla --ref tests/test1.csv  --save {tmp_path} --plot"
    entrypoint(command.split()[1:])

    command = f"star_allele_comp tests/test3.csv --family hla --save {tmp_path}/ --plot"
    proc = subprocess.run(command, shell=True)
    assert proc.returncode == 0
