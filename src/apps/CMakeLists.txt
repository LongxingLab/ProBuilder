# move the final executable file to the compiling folder
set(EXECUTABLE_OUTPUT_PATH ../../)

set( EXES "wefold" "test" "dock" "hallucination" "wefold_ligand_dock" "metalloprotein_sampler" "ligand_binder_sampler" "ligand_binder_sampler_C2" "diversify_peptide" "fold_on_target" "wefold_interface_score" "pH_responsive_assembly" "covid19_vaccine" "interface_metal_heterodimer" "interface_metal_homooligomer" "conformational_switch" "wefold_NCAA" "wefold_NCAA_symm" "GFP_sampler" "light_oligomer" "light_heterodimer" "wefold_pdb_extension" "wefold_motif_pair" "SMIF_homooligomer")
foreach( EXE ${EXES} )

    message( "WeFold exe: " ${EXE} )

    add_executable( ${EXE} ${EXE}.cc args.cc )
    target_link_libraries ( ${EXE} sampling scene utils scoring metric gzip dssp basic)

endforeach( EXE )
