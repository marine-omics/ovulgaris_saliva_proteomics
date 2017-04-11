
# Create bash alias or symbolic link inside new bedtools dir = shortcut to type multiple commands

# Use ln - s to create a symbolic link

```bash
 ln -s ../proteogenomics/novel_pg.gff3 .
```

# use ls -l to list files which will show to which file they are linked to (e.g known.gff3 -> ../proteogenomics/known.gff3)

#bash alias to view colour (so when ls is used it now incorporates "-G")

```bash
alias ls='ls -G'
```
