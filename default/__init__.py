conf_default = {
        "circular": {
            "tickSize": 5
        },
        "features": {
            "fallbackStyle": {
                "color": "#787878",
                "form": "rect",
                "height": 30,
                "visible": False
            },
            "showAllFeatures": False,
            "supportedFeatures": {}
        },
        "graphicalParameters": {
            "buttonWidth": 90,
            "canvasHeight": 900,
            "canvasWidth": 2000,
            "fade": 0.8,
            "genomeLabelWidth": 300,
            "karyoDistance": 5000,
            "karyoHeight": 30,
            "linkKaryoDistance": 20,
            "tickDistance": 1000,
            "tickLabelFrequency": 10,
            "treeWidth": 200
        },
        "labels": {
            "chromosome": {
                "showChromosomeLabels": False
            },
            "features": {
                "showFeatureLabels": False
            },
            "genome": {
                "color": "#000000",
                "showGenomeLabels": True,
                "size": 25
            },
            "showAllLabels": False,
            "ticks": {
                "color": "#000000",
                "showTickLabels": True,
                "showTicks": True,
                "size": 10
            }
        },
        "layout": "linear",
        "linear": {
            "drawAllLinks": False,
            "endLineColor": "#1d91c0",
            "hideHalfVisibleLinks": False,
            "startLineColor": "#1d91c0"
        },
        "maxLinkIdentity": 100,
        "maxLinkIdentityColor": "#D21414",
        "maxLinkLength": 5000,
        "midLinkIdentity": 85,
        "midLinkIdentityColor": "#FFEE05",
        "minLinkIdentity": 60,
        "minLinkIdentityColor": "#1DAD0A",
        "minLinkLength": 100,
        "offset": {
            "distance": 1000,
            "isSet": False
        },
        "tree": {
            "drawTree": True,
            "orientation": "left"
        }
    }
filters_default = {
    "features": {
        "invisibleFeatures": {}
    },
    "karyo": {
        # "chromosomes": {
        #     "seq0": {
        #         "reverse": False,
        #         "visible": True
        #     },
    # },
    "genome_order": [],
    "order": []
},
    "links": {
        "invisibleLinks": {},
        "maxLinkIdentity": 100,
        "maxLinkLength": 301000,
        "minLinkIdentity": 60,
        "minLinkLength": 1000
    },
    "onlyShowAdjacentLinks": True,
    "showAllChromosomes": False,
    "showIntraGenomeLinks": False,
    "skipChromosomesWithoutLinks": False,
    "skipChromosomesWithoutVisibleLinks": True,
}
data_default = {"features":{"link":{}},
                "karyo":{},
                "links":{},
                "tree":{}}

json_obj_default = {'conf':conf_default,
                    'filters':filters_default,
                    'data':data_default}