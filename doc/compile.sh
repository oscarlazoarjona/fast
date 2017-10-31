# A script to compile the documentation.

# We re-install fast to reflect any changes in the source.
pip install .. --upgrade
# We compile the documentation.
make html
