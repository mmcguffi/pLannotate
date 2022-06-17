from plannotate import resources
from plannotate import annotate
from plannotate import streamlit_app

with open("./tests/test_data/RRNB_fragment.txt") as f:
    RRNB = f.read()


def test_yaml():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    assert isinstance(snapgene_db, dict)


def test_BLAST():
    db_meta = resources.get_yaml(resources.get_yaml_path())
    snapgene_db = db_meta["snapgene"]
    hits = annotate.BLAST(RRNB, snapgene_db)
    assert not hits.empty

def test_get_fasta():
    print(resources.get_example_fastas())

def test_databases_exist():
    assert resources.databases_exist() == True
    
def test_valid_sequence_correct():
    resources.validate_sequence("ACTG")
    
def test_valid_sequence_incorrect_base():
    try:
        resources.validate_sequence("ACTX")
        assert False
    except ValueError:
        pass

def test_get_name_ext():
    name, ext = resources.get_name_ext("./a/long/path/test.fasta")
    assert name == "test"
    assert ext == ".fasta"

def test_streamlit_app():
    """this component is hard to test"""
    streamlit_app.run_streamlit(['--yaml-file',resources.get_yaml_path()])
    

def test_annotate():
    hits = annotate.annotate(RRNB)
    assert hits.iloc[0]["sseqid"] == "rrnB_T1_terminator"
    
    