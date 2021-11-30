FROM shinaoka/vscode_julia_python

RUN pip3 install irbasis3 xprec matplotlib
RUN julia -e "import Pkg; Pkg.precompile();"
