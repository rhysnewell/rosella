#!/bin/bash -e

echo "Building Markdown versions of man pages .."
for SUBCOMMAND in recover refine
do
    echo "Converting $SUBCOMMAND .."
    sed 's/\\\[/[/g; s/\\\]/]/g' docs/usage/rosella-$SUBCOMMAND.wd.md |cat <(sed s/SUBCOMMAND/$SUBCOMMAND/ prelude) - >docs/usage/rosella-$SUBCOMMAND.md
    echo "Finished documenting $SUBCOMMAND"

done
rm docs/usage/rosella-*.wd.*
sed -i 's/---NAME//' docs/usage/rosella-*.md