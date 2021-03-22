import VaporSorptionSystem

print("Please enter the name (with \".xlsx\") of the filled template file. \n"
      "If you don't have a template file, please run the template generator executable.")
path_name = input("> ")
system = VaporSorptionSystem.ExcelHandler.read_template_into_sorption_system(path_name)
VaporSorptionSystem.ExcelHandler.write_results(sorption_system=system, output_path="analyzed_" + path_name)
print("Done! Press any key to exit.")
input("")