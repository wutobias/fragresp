install:
	python setup.py install

clean:
	rm -rf *.egg-info
	rm -rf build
	rm -rf dist

clean.pip:
	pip uninstall fragresp --yes

clean.all:
	rm -rf *.egg-info
	rm -rf build
	rm -rf dist
	pip uninstall fragresp --yes
	find -name "*.pyc" -exec rm -f {} \;

clean.example:
	rm -rf example/frag_dir
	rm -rf example/surr_cap_dir
	rm -f example/frag_db.pickle
	rm -f example/*.log
	rm -rf example/datadump
	rm -rf example/remap*
