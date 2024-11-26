from fastapi import FastAPI, UploadFile, File, HTTPException, Query
from fastapi.responses import JSONResponse, RedirectResponse
from typing import Literal
import io
from pydantic import BaseModel
from plannotate.annotate import annotate
from Bio.SeqIO import read
from plannotate import resources as rsc
from plannotate import __version__ as version
import os

MAX_SEQUENCE_LENGTH = (
    os.environ["MAX_SEQUENCE_LENGTH"] if "MAX_SEQUENCE_LENGTH" in os.environ else 20000
)


class ReportRow(BaseModel):
    sseqid: str
    start_location: int
    end_location: int
    strand: int
    percent_identity: float
    full_length_of_feature_in_db: int
    length_of_found_feature: int
    percent_match_length: float
    fragment: bool
    database: str
    Feature: str
    Type: str
    Description: str
    sequence: str


class AnnotateResponse(BaseModel):
    gb_file: str
    report: list[ReportRow]
    version: Literal[f"{version}"]


app = FastAPI()


@app.get("/")
async def root():
    return RedirectResponse(url="/docs")


@app.post("/annotate", response_model=AnnotateResponse)
async def analyze_sequence(
    file: UploadFile = File(...), circular: bool = Query(default=False)
):
    # Read uploaded file contents
    contents = await file.read()
    try:
        if file.filename.endswith((".gbk", ".gb")):
            file_streamer = io.StringIO(contents.decode())
            record = read(file_streamer, "genbank")
        elif file.filename.endswith((".fasta", ".fa")):
            file_streamer = io.StringIO(contents.decode())
            record = read(file_streamer, "fasta")
        elif file.filename.endswith(".dna"):
            file_streamer = io.BytesIO(contents)
            record = read(file_streamer, "snapgene")
        else:
            raise HTTPException(
                status_code=400,
                detail="Invalid file type, only .gbk, .gb, .fasta, .fa, and .dna are supported",
            )
    except Exception:
        raise HTTPException(status_code=400, detail="Error reading file")

    circular = circular or (
        "topology" in record.annotations.keys()
        and record.annotations["topology"] == "circular"
    )
    inSeq = str(record.seq)
    if len(inSeq) > MAX_SEQUENCE_LENGTH:
        raise HTTPException(
            400,
            f"Sequence length is greater than the maximum allowed length of {MAX_SEQUENCE_LENGTH} characters",
        )
    recordDf = annotate(inSeq, is_detailed=True, linear=not circular)
    gbk = rsc.get_gbk(recordDf, inSeq, is_linear=not circular)
    csv = rsc.get_clean_csv_df(recordDf)

    fixed_report = [
        {k.replace(" ", "_"): v for k, v in r.items()}
        for r in csv.to_dict(orient="records")
    ]
    response = {
        "version": version,
        "gb_file": gbk,
        "report": fixed_report,
    }

    return JSONResponse(content=response)
