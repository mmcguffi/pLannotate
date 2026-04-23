import pytest

from plannotate import resources

pytestmark = pytest.mark.integration


def test_streamlit_app():
    """this component is hard to test"""
    from plannotate import streamlit_app

    streamlit_app.run_streamlit(["--yaml-file", resources.get_yaml_path()])


# # runs indefinitely
# def test_streamlit():
#     runner = CliRunner()
#     result = runner.invoke(main_streamlit)
#     assert result.exit_code == 0
