"""
Unit tests for EuropePMCClient and snippet extraction.
"""

from unittest.mock import Mock, patch
import pytest

from magnamined.literature import EuropePMCClient


def test_extract_targeted_snippets():
    text = (
        "Background: Rice blast is caused by Magnaporthe oryzae. " * 10
        + "Results: The MGG_01001 deletion mutant exhibited severe virulence reduction on rice leaves. "
        + "Complementation rescued the wild-type phenotype. "
        + "Discussion: This demonstrates essential role in infection. " * 10
    )

    snippets = EuropePMCClient.extract_targeted_snippets(
        text_content=text,
        gene_id="MGG_01001",
        aliases=["ENG1"],
        window_chars=300
    )

    assert len(snippets) > 0
    assert any("MGG_01001" in s for s in snippets)
    assert any("virulence reduction" in s for s in snippets)


def test_search_gene_query_construction():
    client = EuropePMCClient()
    with patch.object(client.session, "get") as mock_get:
        mock_resp = Mock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {"resultList": {"result": []}}
        mock_get.return_value = mock_resp

        client.search_gene("MGG_01001", gene_name="ENG1")

        args, kwargs = mock_get.call_args
        query_param = kwargs["params"]["query"]
        assert '"MGG_01001"' in query_param
        assert '"ENG1"' in query_param
        assert '"MoENG1"' in query_param
        assert 'Magnaporthe oryzae' in query_param
