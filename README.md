# FPC_mRx

## Overview
This repository contains a complete reimplementation of the Magnetic Reconnection analysis pipeline originally developed in IDL by Arya. The new version is written in Python for modularity, transparency, and long-term maintainability.

## Structure
- `/fpc_mrx_idl_archive`: Legacy IDL scripts and original documentation
- `/src`: Core Python modules for data loading, processing, and analysis
- `/notebooks`: Exploratory notebooks for validation, visualization, and prototyping
- `/data`: Sample input/output formats and test datasets

## Transition Notes
The IDL scripts have been archived for reference only in the folder `fpc_mrx_idl_archive`. All logic has been restructured in Python with updated conventions, modular design, and improved readability. No direct compatibility is maintained.

## Installation
Clone the repository and install dependencies:
```bash
git clone git@research-git.uiowa.edu:howes-group/fpc_mrx.git
cd fpc_mrx
pip install -r requirements.txt
```

## Author
Regis John  
Postdoctoral Researcher, University of Iowa  
Specializing in magnetic reconnection and plasma dynamics
