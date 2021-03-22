import PyInstaller.__main__

generator_exe_name = 'generate_vapor_sorption_template'
generator_script_name = 'generate_template.py'

analyzer_exe_name = 'analyze_vapor_sorption_template'
analyzer_script_name = 'analyze_template.py'


print("Compiling template generator...")
PyInstaller.__main__.run([
    '--name=%s' % generator_exe_name,
    '--onefile',
    '--windowed',
    '--hidden-import=%s' % 'pkg_resources.py2_warn',  # typically PyInstaller fails to find this and gives a loading error, so we will force it to load here.
    #'--icon=%s' % icon_path,
    generator_script_name,
])
print("Finished compiling template generator.")


print("Compiling template analyzer...")
PyInstaller.__main__.run([
    '--name=%s' % analyzer_exe_name,
    '--onefile',
    # '--windowed',
    '--console',
    '--hidden-import=%s' % 'pkg_resources.py2_warn',  # typically PyInstaller fails to find this and gives a loading error, so we will force it to load here.
    #'--icon=%s' % icon_path,
    analyzer_script_name,
])
print("Finished compiling template analyzer.")

