def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num - 1] = text + '\n'  # array index starts at 0, subtract 1
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()


first = 52
last  = 96
for inputNumber in range(first,last+1):
    replace_line('MAIN.py', 20, f"exec(open('Input/inputData{inputNumber}.py').read())")
    replace_line('MAIN.py', 95, f"numFolder = {inputNumber}")
    exec(open("MAIN.py").read())






















