# A script to compile the documentation.

# We re-install fast to reflect any changes in the source.
pip install .. --upgrade
# We go to the documentation repository:
cd ../../website/html
# We empty it.
git clean -fdx

# We go back to the fast reposotory
cd ../../fast/doc

# We compile the documentation.
make html

# We go to the documentation repository:
cd ../../website/html
# We add the changes.
git add .
# Commit
git commit -m "Edit."
# We can push at this point if the result is good with.
# git push origin gh-pages
