import subprocess
from pathlib import Path


def test_blastn_builds_command_and_parses_tabular_output(monkeypatch):
    from plannotate import annotate

    commands = []

    def fake_run(args, shell, capture_output, text):
        commands.append(args)
        output_path = Path(args[args.index("-out") + 1])
        output_path.write_text("1 4 feature 1 100 4 ACTG 4 1 4 4 0.0\n")
        return subprocess.CompletedProcess(args, 0)

    monkeypatch.setattr(annotate.subprocess, "run", fake_run)

    hits = annotate.BLAST(
        "ACTG",
        {
            "method": "blastn",
            "parameters": "-perc_identity 95",
            "db_loc": "/db/snapgene",
        },
    )

    command = commands[0]
    assert command[0] == "blastn"
    assert command[1:3] == ["-task", "blastn-short"]
    assert command[command.index("-db") + 1] == "/db/snapgene"
    assert "-perc_identity" in command
    assert "6 qstart qend sseqid" in command[command.index("-outfmt") + 1]
    assert hits.loc[0, "sseqid"] == "feature"
    assert hits.loc[0, "qseq"] == "ACTG"
    assert hits.loc[0, "evalue"] == 0.0


def test_diamond_builds_command_and_normalizes_subject_ids(monkeypatch):
    from plannotate import annotate

    commands = []

    def fake_run(args, shell, capture_output, text):
        commands.append(args)
        output_path = Path(args[args.index("-o") + 1])
        output_path.write_text("1 4 sp|PROT|desc 100 2 MAAA 4 1 2 4 0.0\n")
        return subprocess.CompletedProcess(args, 0)

    monkeypatch.setattr(annotate.subprocess, "run", fake_run)

    hits = annotate.BLAST(
        "ATGGCGGCGGCG",
        {
            "method": "diamond",
            "parameters": "--id 75",
            "db_loc": "/db/fpbase",
        },
    )

    command = commands[0]
    assert command[:2] == ["diamond", "blastx"]
    assert command[command.index("-d") + 1] == "/db/fpbase"
    assert "--id" in command
    outfmt_index = command.index("--outfmt")
    assert command[outfmt_index + 1 : outfmt_index + 5] == [
        "6",
        "qstart",
        "qend",
        "sseqid",
    ]
    assert hits.loc[0, "sseqid"] == "PROT"
    assert hits.loc[0, "sframe"] == 1
    assert hits.loc[0, "slen"] == 6
    assert hits.loc[0, "length"] == 4
