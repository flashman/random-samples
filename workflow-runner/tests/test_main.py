import asyncio
import re
import pytest
from workflow_runner.main import WorkflowRunner, DAG


@pytest.fixture
def json_input_super_simple():
    return """
    {
        "A": {"start": true, "edges": {"B": 1}},
        "B": {"edges": {}}
    }
    """


@pytest.fixture
def json_input_simple():
    return """
    {
        "A": {"start": true, "edges": {"B": 2, "C": 1}},
        "B": {"edges": {"D": 3}},
        "C": {"edges": {"D": 5}},
        "D": {"edges": {}}
    }
    """


@pytest.fixture
def json_input_multi_visit():
    return """
    {
        "A": {"start": true, "edges": {"B": 2, "C": 1}},
        "B": {"edges": {"D": 3}},
        "C": {"edges": {"D": 5}},
        "D": {"edges": {"E": 1}},
        "E": {"edges": {}}
    }
    """


@pytest.fixture
def json_input_slow_level():
    return """
    {
        "A": {"start": true, "edges": {"B": 1, "C": 5}},
        "B": {"edges": {"D": 2}},
        "C": {"edges": {}},
        "D": {"edges": {"E": 1}},
        "E": {"edges": {}}
    }
    """


@pytest.fixture
def json_input_cycle():
    return """
    {
        "A": {"start": true, "edges": {"B": 2}},
        "B": {"edges": {"C": 1}},
        "C": {"edges": {"A": 1}}
    }
    """


@pytest.fixture
def json_input_no_start():
    return """
    {
        "A": {"edges": {"B": 2}},
        "B": {"edges": {"C": 1}},
        "C": {"edges": {"D": 3}},
        "D": {"edges": {}}
    }
    """


@pytest.fixture
def json_input_missing_node():
    return """
    {
        "A": {"edges": {"B": 1}},
    }
    """


def test_DAG_with_cycle(json_input_cycle):
    with pytest.raises(ValueError):
        DAG.model_validate_json(json_input_cycle)


def test_DAG_no_start(json_input_no_start):
    with pytest.raises(ValueError):
        DAG.model_validate_json(json_input_no_start)


def test_DAG_missing_node(json_input_missing_node):
    with pytest.raises(ValueError):
        DAG.model_validate_json(json_input_missing_node)


@pytest.mark.asyncio
async def test_workflow_runner_super_simple(json_input_super_simple, capsys):
    runner = WorkflowRunner(json_input_super_simple)
    await runner.run()

    captured = capsys.readouterr()
    assert captured.out.replace("\n", "") == "AB"


@pytest.mark.asyncio
async def test_workflow_runner_simple(json_input_simple, capsys):
    runner = WorkflowRunner(json_input_simple)
    await runner.run()

    captured = capsys.readouterr()
    assert captured.out.replace("\n", "") == "ACBD"


@pytest.mark.asyncio
async def test_workflow_runner_multi_visit(json_input_multi_visit, capsys):
    runner = WorkflowRunner(json_input_multi_visit)
    await runner.run()

    captured = capsys.readouterr()
    assert captured.out.replace("\n", "") == "ACBDE"


@pytest.mark.asyncio
async def test_workflow_runner_slow_level(json_input_slow_level, capsys):
    runner = WorkflowRunner(json_input_slow_level)
    await runner.run()

    captured = capsys.readouterr()
    assert captured.out.replace("\n", "") == "ABDEC"
