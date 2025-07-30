import numpy as np
import pytest
import sys
from pathlib import Path

from julesf.coupler.coupler import run_coupled_model
from julesf.rothc.parameters import POOLS, C0_default

def test_coupled_model_basic():
    triffid_init = [6.0, 3.0, 0.8, 0.2]
    rothc_init = [C0_default[p] for p in POOLS]
    
    initial_conditions = {
        'triffid': triffid_init,
        'rothc': rothc_init
    }
    
    drivers_ts = {
        'T_soil': lambda t: 283.15,
        'moisture': lambda t: 0.5,
        's': lambda t: 0.5
    }
    
    results = run_coupled_model(
        t_span_weeks=(0, 10), 
        initial_conditions=initial_conditions
    )
    
    assert 'time_weeks' in results   
    assert 'triffid' in results
    assert 'rothc' in results
    assert results['triffid'].shape[1] == len(results['time_weeks'])
