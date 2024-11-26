from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.responses import JSONResponse, RedirectResponse
import io
from pydantic import BaseModel
from plannotate.annotate import annotate
from Bio.SeqIO import read
from plannotate import resources as rsc


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


app = FastAPI()


@app.get("/")
async def root():
    return RedirectResponse(url="/docs")


@app.post("/annotate", response_model=AnnotateResponse)
async def analyze_sequence(file: UploadFile = File(...)):
    # Read uploaded file contents
    contents = await file.read()

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
        raise HTTPException(status_code=400, detail="Invalid file type")

    inSeq = str(record.seq)
    recordDf = annotate(inSeq, is_detailed=True, linear=False)
    gbk = rsc.get_gbk(recordDf, inSeq, is_linear=False)
    csv = rsc.get_clean_csv_df(recordDf)
    print(gbk)
    response = {
        "gb_file": gbk,
        "report": csv.to_dict(orient="records"),
    }

    return JSONResponse(content=response)
