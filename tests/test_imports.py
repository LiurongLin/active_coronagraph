def test_package_imports():
    import active_coronagraph
    import active_coronagraph.basic
    import active_coronagraph.cli
    import active_coronagraph.coronagraphs
    import active_coronagraph.coronagraphs_polished
    import active_coronagraph.main_functions_polished
    import active_coronagraph.make_plot
    import active_coronagraph.new_mask
    import active_coronagraph.phase_masks

    assert active_coronagraph.__version__
