project('Pr�ctica 2 de PAV - detecci�n de actividad vocal', 'c')

exe = 'vad'
src = ['src/main_vad.c', 'src/vad.c', 'src/pav_analysis.c']
lib = ['-lm', '-lsndfile', '-lgsl', '-lgslcblas']

executable(exe, sources: src, link_args: lib)
    