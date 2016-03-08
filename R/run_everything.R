## Copyright (c) 2016, Otávio Bartalotti, Gray Calhoun, and Yang He.
## Available under the MIT "Expat" License, see README.md

## Like the file name says, this file runs everything. Each file is
## run in its own environment so we avoid clutter, etc.

sys.source("replication.R", environment())
sys.source("application.R", environment())
sys.source("simulation.R", environment())
