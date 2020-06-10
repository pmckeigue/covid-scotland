# covid-scotland


## Getting started with Github, from a Linux command line


### Saving your credentials in a config file (~/.git-credentials)

> `git config --global credential.helper store"`

The first time you pull files from the repository you will be asked for a username and password, but after that they will be stored and you won't have to type them again. 

### Make a local copy of the repository

> `git clone https://github.com/pmckeigue/covid-scotland.git`


### Add a local file to the list of files under Github version control

> `git add myscript.R`

### Commit the latest version of this file before sending it

> `git commit -m"your message about this update" myscript.R`


### Send all latest commits to the repository

> `git push -u origin master`

(the -u option adds a tracking reference that can be used by git pull)

### View status of all local files

> `git status`

### Update the local copy from the repository

> `git pull origin master`

This is a useful tutorial on using git from the command line

https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html#add-and-commit-local-changes
