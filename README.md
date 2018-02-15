# Cluster viewer
Visualize cluster data from ChIB or WFA experiments in browser.
Data should be in bed format with cluster id in the forth column. Example:
```
chr16	71461594	71461783	33254	2	+
chr2	210093832	210093994	33254	8	+
chr7	128739812	128740098	33010	10	+
...
```

To start bokeh server move into directory:
```
cd /path/to/Cluster_viewer
```

and type the following command
``` 
bokeh serve . --show
```

Data is required to be located or linked-to in the "data" folder.
Alternatively the following command can be used to input data directly.
```
bokeh serve . --show --arg /path/to/bed_file
```

Copy the URL e.g (http://localhost:5006/Cluster_viewer or similar) into your browser.

Good to go!!

