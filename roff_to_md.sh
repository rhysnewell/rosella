#!/bin/bash -e

echo "Building Markdown versions of man pages .."
for SUBCOMMAND in recover refine
do
    echo "Converting $SUBCOMMAND .."
    # Demote pandoc's man sections to H2 so the prelude heading is the page's
    # only H1, then drop the now-redundant NAME section header.
    sed -e 's/\\\[/[/g; s/\\\]/]/g' \
        -e 's/^# /## /' \
        -e '/^## NAME$/d' \
        docs/usage/rosella-$SUBCOMMAND.wd.md \
      | cat <(sed "s/SUBCOMMAND/$SUBCOMMAND/" prelude) - \
      > docs/usage/rosella-$SUBCOMMAND.md
    echo "Finished documenting $SUBCOMMAND"
done
rm docs/usage/rosella-*.wd.*
