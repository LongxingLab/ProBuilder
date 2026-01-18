# move the final executable file to the compiling folder
set(EXECUTABLE_OUTPUT_PATH ../../)

set( EXES "inducible_homooligomer_generator" "scaffold_generator" "light_heterodimer" "light_homooligomer" "covalent_ligand_binder")
foreach( EXE ${EXES} )

    message( "WeFold exe: " ${EXE} )

    add_executable( ${EXE} ${EXE}.cc args.cc )
    target_link_libraries ( ${EXE} sampling scene utils scoring metric gzip dssp basic)

endforeach( EXE )
