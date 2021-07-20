if (FALSE) {
  igv <- file.path("c:", "Program Files", "IGV_2.8.13", "igv.bat")

  shell("java --module-path=lib -Xmx4g @igv.args --module=org.igv/org.broad.igv.ui.Main")

  shell("java --module-path=lib -Xmx4g --module=org.igv/org.broad.igv.ui.Main")

  shell.exec(igv)

  bamfolder <-
    file.path("c:", "Users", "steve rozen", "Desktop", "bamSlice_folder")

  Tbam.name <-
    dir(file.path(bamfolder, "DO52605_SP116496_dbs_srt_fromBED.bam"),
        full.names = TRUE)

  Nbam.name <-
    dir(file.path(bamfolder, "DO52605_SP116498_dbs_srt_fromBED.bam"),
        full.names = TRUE)

}


create_IGV_snapshot_script <-
  function(vcf.name,
           Tbam.name,
           Nbam.name,
           igv.script.name = paste0(vcf.name, "_igv_script.txt"),
           out.dir = paste0(vcf.name, "_igv_script_out"),
           genome = "Human hg19") {

    if (!dir.exists(out.dir)) {
      if (!dir.create(out.dir)) {
        warning("Did not create", out.dir)
      }
    }

    pp <- function(...) {
      cat(..., "\n", sep = "", file = igv.script.name, append = TRUE)
    }



    cat("new\n", file = igv.script.name)
    pp("genome ", genome)
    pp("load ", Tbam.name)
    pp("Load ", Nbam.name)
    pp("snapshotDirectory ", out.dir)



  }
