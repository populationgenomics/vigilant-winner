# PanelApp

PanelApp is a crowd-sourced gene curation knowledge base, aggregating information on gene-disease associations from a wide range of specialists. Within PanelApp there are a number of panels, each containing a set of genes associated with a disease, or a theme. Within a panel genes are rated using a traffic light system:

* Red (minimal evidence)
* Orange (intermediate evidence)
* Green (strong evidence linking gene to theme)

For this project, we will focus analysis on a specific list of genes associated strongly with mendelian disease, with an option to supplement this gene list with any number of additional PanelApp entries. The Gene list we will use as the core of the analysis defaults to the [Mendeliome](https://panelapp-aus.org//panels/137/), though this can be altered through the configuration file.

During an analysis run, the [Query Script](../src/talos/QueryPanelapp.py) is used to pull down the latest core panel data (gene ENSG, symbol, and Mode of Inheritance), saved as a dictionary indexed on GRCh38/Ensembl 90 ENSG.

## Per-participant Panels

This codebase contains a [HPO-Panel Helper](../src/talos/GeneratePanelData.py), integrating participant phenotypic data, PanelApp per-panel relevant phenotypes, and the HPO heirarchy to find relevant panels to apply to each individual participant. This will be difficult to replicate outside of CPG infrastructure, so this describes the output generated, which can be produced by other means:

```python
# participantID here must match the ID in the joint call
_panel_dict = {
    'participantID': {
        'panels': [1, 23, 456]
    }
}
```

If this file is supplied at runtime, any panel matched to a participant will be integrated into the overall panel query. Algorithmically, the following process occurs:

1 Query for the core panel (Mendeliome, unless altered in config)

* We collect the relevant MOI and symbol for each Green gene

2 Gather all phenotype-matched panels from the file
3 Query for each extra panel in turn, adding new genes from the additional panels

* At each stage, if MOI was missing in a prior panel and present in the current panel, update the MOI
* Record all the panels each gene appears on
