# Minimal makefile for Sphinx documentation
# e.g. make html

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = docsrc
BUILDDIR      = docbuild


# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

doc:
	@make html
	@rm -rf ./docs/*
	@mkdir -p ./docs
	@touch ./docs/.nojekyll
	@cp -a docbuild/html/. ./docs
	@git add docs/
	@rm -rf docbuild

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)


# Building python
clean:
	@rm -rf *.egg-info/
	@rm -rf dist/
	@rm -rf build/

tests:
	@python3 -m unittest discover tests/

upload:
	@python3 -m twine upload dist/* --verbose

upload-dev:
	@python3 -m twine upload --repository testpypi --verbose dist/*

build:
	@git branch --color=never | grep "* master" > /dev/null
	@make tests
	@make clean
	@make doc
	@python3 setup.py sdist bdist_wheel
	@rm -rf *.egg-info/
	@rm -rf build/

release:
	@make build
	@make upload

release-dev:
	@make build
	@make upload-dev
