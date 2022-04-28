#' @export
gimmefMRI <- function(
      datadir = 'C:/Users/j186b025/Documents/GitHub/jbartolotti/gimmefMRI/demodat',
      savedir = '//kumc.edu/data/Research/Hoglund/Bartolotti_J/gimme_toolkit/models',
      run_now = TRUE
      ){

  mm <- readModelInput(file.path(datadir,'timecourses.csv'),
                 file.path(datadir,'model_spec.csv'),
                 file.path(datadir,'info_lists.csv'),
                 file.path(datadir,'shortnames.csv'))

  #Creates model folders and fills input_files folders with
  #subject timecourses from mm$timecourses based on
  #model parameters in mm$model_spec
  initializeGimmeFolders(savedir, mm)

  #Creates an .R file containing code to run all models in parallel
  writeGimmeCode(savedir, mm)
#  if(run_now){
#    runGimme()
#  }
}


