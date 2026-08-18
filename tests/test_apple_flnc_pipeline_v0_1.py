import gzip
import json
import sys
from pathlib import Path
from unittest.mock import patch


SCRIPT = (
    Path(__file__).parents[1]
    / "scripts"
    / "run_apple_golden_delicious_ont_flnc_validation_v0.1.sh"
)


def embedded_python_programs() -> list[str]:
    chunks = SCRIPT.read_text(encoding="utf-8").split("<<'PY'\n")[1:]
    programs = []
    for chunk in chunks:
        source, separator, _ = chunk.partition("\nPY\n")
        assert separator
        programs.append(source)
    return programs


def test_flnc_pipeline_freezes_published_protocol_and_strand_semantics() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "PLOIDYPATCH_CODE_COMMIT must be a full git SHA" in text
    assert "[[ $pychopper_version == 2.7.10 ]]" in text
    assert "--minimum-length 500" in text
    assert "--minimum-length 500 --minimum-mean-quality" not in text
    assert ".len500.fastq.gz" in text
    assert ".q7_len500.fastq.gz" not in text
    assert "mean_quality_filter\\tPychopper_mean_error_probability_PHRED_ge_7" in text
    assert "-k PCS109 -m phmm -Q 7 -z 50 -t 16" in text
    assert "autotune_seed=20261005" in text
    assert '[[ ${pandas_version%%.*} == 2 ]]' in text
    assert "Pychopper_2.7.10_reporting_requires_pandas_2.x" in text
    assert "Bioconda_curated_modern_runtime_metadata_exception" in text
    assert "SOURCE_DATE_EPOCH=0" in text
    assert "source_date_epoch\\t0" in text
    assert "def count_fastq_records(path):" in text
    assert '"pychopper_primary_classified": primary_classified' in text
    assert '"pychopper_primary_flnc": primary_flnc' in text
    assert '"pychopper_output_segments_below_50": output_length_failed' in text
    assert '"pychopper_rescued_reads_excluded": rescued_reads' in text
    assert 'PATH="$pychopper_env/bin:/usr/bin:/bin"' in text
    assert "-x splice -uf -k14 -G 1000000" in text
    assert text.count("--alignment-strand-source query_orientation") == 2
    assert "--alignment-strand-source minimap2_ts" not in text
    assert "pychopper_rescued_segments_used\\tfalse" in text


def test_flnc_pipeline_is_non_overwriting_and_retains_negative_results() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "refusing to overwrite apple Pychopper FLNC validation" in text
    assert "candidate_coordinates_modified\\tfalse" in text
    assert "candidate_ranks_modified\\tfalse" in text
    assert "evidence_thresholds_modified\\tfalse" in text
    assert "reported_ranges_are_descriptive_not_validity_gates" in text
    count_gate = "if not 0 < primary_flnc <= primary_classified <= pass_reads"
    assert count_gate in text
    assert "case_study_ready" not in text.split(count_gate, 1)[0]
    assert "raw_ts_patch2_vs_flnc.json" in text
    assert "chmod -R a-w" in text
    assert "mv \"$working_root\" \"$result_root\"" in text


def test_all_embedded_python_programs_compile() -> None:
    programs = embedded_python_programs()
    assert len(programs) == 3
    for index, source in enumerate(programs, start=1):
        compile(source, f"{SCRIPT.name}:heredoc-{index}", "exec")


def test_classification_summary_counts_materialized_flnc_records(tmp_path: Path) -> None:
    contract = tmp_path / "contract.tsv"
    filter_root = tmp_path / "filter"
    flnc_root = tmp_path / "flnc"
    filter_root.mkdir()
    flnc_root.mkdir()
    rows = [f"ACC{index}\ttissue{index}\tmd5\turl" for index in range(7)]
    contract.write_text(
        "accession\ttissue\texpected_md5\turl\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    stats = (
        "Category\tName\tValue\n"
        "ReadStats\tPassReads\t4\n"
        "ReadStats\tQcFail\t0\n"
        "ReadStats\tLenFail\t1\n"
        "Classification\tPrimers_found\t2\n"
        "Classification\tRescue\t2\n"
        "Classification\tUnusable\t1\n"
        "RescueSegmentNr\t2\t1\n"
    )
    for index in range(7):
        accession = f"ACC{index}"
        (filter_root / f"{accession}.filter_summary.json").write_text(
            json.dumps(
                {"input_records": 4, "passed_records": 4, "length_failed_records": 0}
            ),
            encoding="utf-8",
        )
        (flnc_root / f"{accession}.stats.tsv").write_text(stats, encoding="utf-8")
        with gzip.open(flnc_root / f"{accession}.flnc.fastq.gz", "wt") as handle:
            handle.write("@read\nACGT\n+\nIIII\n")
    output_tsv = tmp_path / "summary.tsv"
    output_json = tmp_path / "summary.json"
    argv = [
        "classification-summary",
        str(contract),
        str(filter_root),
        str(flnc_root),
        str(output_tsv),
        str(output_json),
    ]
    with patch.object(sys, "argv", argv):
        namespace = {"__name__": "__main__"}
        exec(
            compile(embedded_python_programs()[0], "classification-summary", "exec"),
            namespace,
            namespace,
        )
    report = json.loads(output_json.read_text(encoding="utf-8"))
    assert report["totals"]["pychopper_primary_classified"] == 14
    assert report["totals"]["pychopper_primary_flnc"] == 7
    assert report["totals"]["pychopper_output_segments_below_50"] == 7
    assert report["totals"]["pychopper_rescued_reads_excluded"] == 7
    assert report["totals"]["pychopper_rescued_segments_excluded"] == 14
