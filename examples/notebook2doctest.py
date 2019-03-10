# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazo@correo.nucleares.unam.mx
r"""A script to transform jupyter notebooks into doctests."""
import os


def find_braces(s, ini, end):
    r"""This function returns a list of pairs of matching text indices in s
    for parenthesis, braces,  brakets, etc as specified by ini and end.
    """
    toret = {}
    pstack = []

    for i, c in enumerate(s):
        if c == ini:
            pstack.append(i)
        elif c == end:
            if len(pstack) == 0:
                raise IndexError("No matching closing parens at: " + str(i))
            toret[pstack.pop()] = i

    if len(pstack) > 0:
        raise IndexError("No matching opening parens at: " + str(pstack.pop()))

    pairs = []
    for ini in toret:
        pairs += [(ini, toret[ini])]
    pairs.sort()

    return pairs


def find_fields(text, left, right, preceded):
    r"""This function returns a list of instances of strings found in
    "text" that are between left and right one-character markers
    "left" and "right", as long as the left  marker is preceded by
    "preceded".
    """
    brackets = find_braces(text, left, right)
    fields = []
    for ini, end in brackets:
        if text[ini-len(preceded):ini] == preceded:
            fields += [text[ini:end+1]]
    return fields


def get_cells(text):
    r"""This function returns a list of tuples of the form
    (cell_type, source, outputs)
    """
    # We find all instances of opening and closing braces.
    braces = find_braces(text, "{", "}")
    cells = []

    # We determine which pairs of braces correspond to cells.
    for i, (ini, end) in enumerate(braces):
        content = text[ini:end+1]
        # If the opening brace is preceded by cell_ini, then that pair
        # corresponds to a cell.
        start = content[2:2+len(cell_ini)]
        if start == cell_ini:
            cell = content
            cell_type = cell[2+len(cell_ini):cell.find('",')]

            # We get the source of the cell.
            source = find_fields(cell, "[", "]", '"source": ')[0]

            # There may or may not be outputs to the cell.
            outputs = find_fields(cell, "[", "]", '"outputs": ')
            text_output = []
            if len(outputs) > 0:
                text_output += find_fields(outputs[0], "[", "]",
                                           '"text/plain": ')
                text_output += find_fields(outputs[0], "[", "]",
                                           '"text": ')

            cells += [(cell_type, source, text_output)]
    return cells


########################################################################
# The path to the notebooks.
path_notebooks = r"../../fast-notebooks/"

# The notebooks to be converted.
aux = r"04 - Vectors in the helicity basis and the electric field.ipynb"
notebooks = [r"01 - Two level atom.ipynb",
             r"02 - Three level atom (ladder).ipynb",
             r"03 - Rb87 one photon.ipynb",
             aux,
             r"05 - Two level atom (symbolic).ipynb",
             r"06 - Three level atom, ladder (symbolic).ipynb",
             r"07 - Three level atom, Lambda (symbolic).ipynb",
             r"08 - Three level atom, V (symbolic).ipynb",
             r"09 - Thermal States.ipynb",
             r"10 - States database.ipynb"]
########################################################################
# We choose a notebook.
index = 3
notebook_name = notebooks[index-1]
f = file(path_notebooks+notebook_name, "r")
text = f.read()
f.close()

cell_ini = '   "cell_type": "'

cells = get_cells(text)
doc = '''# -*- coding: utf-8 -*-
# Copyright (C) 2017 Oscar Gerardo Lazo Arjona
# mailto: oscar.lazo@correo.nucleares.unam.mx

__doc__ = r"""\n'''
Nt = 46
for i, cell in enumerate(cells[:]):
    cell_code = ""

    if cell[0] == "markdown":
        doc += "\n"
        lines = cell[1].split("\n")
        last_line = ""
        for line in lines:
            line = line.replace('\\"', '"')
            if line not in ['[', '   ]']:
                if line[:5] == '    "':
                    line = line[5:]

                if line[-4:] == '\\n",': line = line[:-4]
                if line[-1:] == '"': line = line[:-1]
                line += '\n'

                cell_code += line
                last_line = line

        if last_line[:4] == "... ":
            cell_code += "... \n"

    elif cell[0] == "code":
        # We convert the source code into >>> and ... entries.
        doc += "\n"
        lines = cell[1].split("\n")
        last_line = ""
        for line in lines:
            line = line.replace('\\"', '"')
            if line not in ['[', '   ]']:

                if line[-4:] == '\\n",': line = line[:-4]
                if line[-1:] == '"': line = line[:-1]

                if line[:5] == '    "':
                    if line[5:9] in ["    ", "else"]:
                        line = "... "+line[5:]
                    else:
                        line = ">>> "+line[5:]

                if "%matplotlib inline" in line:
                    line = ""

                if "pyplot.semilogx" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP1"
                if "pyplot.plot" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP1"
                if "ax.plot" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP1"

                if "pyplot.ylabel" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP2"
                if "pyplot.xlabel" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP2"
                if "pyplot.legend" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP2"
                if "ax.text(" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP2"
                if "ax.legend" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP2"

                if "pyplot.axis" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP3"
                if "ax.set_xlim" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP3"
                if "ax.set_ylim" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP3"
                if "pyplot.xlim" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP3"
                if "pyplot.ylim" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP3"

                if "pyplot.savefig" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP4"
                if "fancy_matrix_plot(" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP4"
                if "fancy_r_plot(" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP4"
                if "draw_lasers_3d(" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP4"
                if "excitation(" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP4"
                if "decay(" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP4"

                if "ax.set_xlabel" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP5"
                if "ax.set_ylabel" in line:
                    line = line+" # doctest: +IGNORE_PLOT_STEP5"

                line = line.replace(r"\\", "\\")
                line += "\n"

                cell_code += line
                last_line = line

        if last_line[:4] == "... ":
            cell_code += "... \n"

        # We enter the results, if any.
        outputs = cell[2]

        for output in outputs:
            # print output
            lines = output.split("\n")
            for line in lines:
                if line not in ['[', '   ]']:
                    line = line.replace(r"\\", "\\")
                    if line[-4:] == '\\n",': line = line[:-4]
                    if line[-1:] == '"': line = line[:-1]
                    if line[:7] == '      "': line = line[7:]
                    if line[:8] == '       "': line = line[8:]
                    if line[-2:] == '\\n': line = line[:-2]
                    if line[-7:] == '      ]' and len(line) == len('      ]'):
                        line = "\n"
                    if line[-6:] == '     ]' and len(line) == len('     ]'):
                        line = "\n"

                    blankline = True
                    for char in line:
                        if char != ' ':
                            blankline = False
                            break
                    if blankline: line = "<BLANKLINE>"

                    line += "\n"
                    cell_code += line

    doc += cell_code

doc += '\n"""\n'
aux = '''__doc__=__doc__.replace("+IGNORE_PLOT_STEP'''
doc += aux+'''1", "+ELLIPSIS\\n[<...>]")\n'''
doc += aux+'''2", "+ELLIPSIS\\n<...>")\n'''
doc += aux+'''3", "+ELLIPSIS\\n(...)")\n'''
doc += aux+'''5", "+ELLIPSIS\\nText(...)")\n'''
doc += aux+'''4", "\\n")\n'''


doctest_name = "doctest_"+notebook_name
doctest_name = doctest_name.replace(' ', '_')
doctest_name = doctest_name.replace('-', '_')
doctest_name = doctest_name.replace(".ipynb", ".py")
doctest_name = doctest_name.replace("(", "")
doctest_name = doctest_name.replace(")", "")
doctest_name = doctest_name.replace(",", "")

f = file(doctest_name, "w")
f.write(doc)
f.close()

try:
    os.system("mkdir "+doctest_name[:-3].replace("doctest", "folder"))
except:
    pass

exec("from doctest import testmod")
s = "import "+doctest_name[:-3]
exec(s)
s = "print testmod("+doctest_name[:-3]+")"
print s
exec(s)


# Rules:
#       cells that return text must do so at the end.
#       cells for plotting must end with a savefig call
#       lines of code cannot end in "
#       pyplot. lines cannot be commented.
