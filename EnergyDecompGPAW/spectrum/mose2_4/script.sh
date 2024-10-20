#!/bin/sh
python3 fdm.py >& fdm.log
python3 ksd.py >& ksd.log
python3 tcm.py >& tcm.log
python3 indden.py >& indden.log
