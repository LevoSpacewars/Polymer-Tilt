import os
from subprocess import call
files = os.listdir(".")

cpp = []
c = []

for i in range(len(files)):
    if(files[i].endswith(".cpp")):
        cpp.append(files[i].strip(".cpp"))
    if(files[i].endswith(".c")):
        c.append(files[i].strip(".c"))

gcc_compile_string = ["gcc","-c","-std=c11","-o","",""]
for i in range(len(c)):
    gcc_compile_string[-2] = c[i] + ".c"
    gcc_compile_string[-1] =  c[i] + ".o"
    call(gcc_compile_string)

print(gcc_compile_string)

gpp_compile_string = ["g++-8", "-c", "-lstdc++fs","-std=c++17", "-o","","" ,"-lstdc++fs"]
for i in range(len(cpp)):
    gpp_compile_string[-3] = cpp[i] + ".o"
    gpp_compile_string[-2] = cpp[i] + ".cpp"
    call(gpp_compile_string)

print(gpp_compile_string)


gpp_linker_string = ["g++-8","-o", "Executable"]
for i in range(len(cpp)):
    gpp_linker_string.append(cpp[i] + ".o")

for i in range(len(c)):
    gpp_linker_string.append(c[i] + ".o")
gpp_linker_string.append("-lstdc++fs")
call(gpp_linker_string)




print(gpp_linker_string)
