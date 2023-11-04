#!/bin/bash -e

echo "Building ROFF versions of man pages .."
for SUBCOMMAND in recover refine
do
    echo "Documenting $SUBCOMMAND .."
    cargo run -- $SUBCOMMAND --full-help-roff > docs/usage/rosella-$SUBCOMMAND.wd.roff
    echo "Finished documenting $SUBCOMMAND"
done