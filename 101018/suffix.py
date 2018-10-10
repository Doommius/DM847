from suffix_trees import STree
st = STree.STree("teethermometer")
print(st.find("mom")) # 0
print(st.find_all("et")) # [0, 8]



string = "teethermometer$"
print (string)





while string != "":
    string = (string[1:])
    print(string)