#!/bin/bash

export SPHINX_APIDOC_OPTIONS='members,show-inheritance,private-members,special-members'
sphinx-apidoc -efMT -o . ../stability_conditions