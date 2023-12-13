#!/bin/bash
g++ -O3 -std=c++20 -DNDEBUG clauses.cpp restore_list.cpp main.cpp -o dynex-testnet-bnb
strip dynex-testnet-bnb 2>/dev/null
