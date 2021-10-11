#Clone the repository to your local computer as follows (substitute username for your username)
git clone username@idocdr1:/storage/cdr/Shared/BioFig.git BioFig.git

#check URLs of remote
git remote -v

#Check what changed in the last commits
git whatchanged origin/master -n 1

git fetch origin
git diff origin/master
git merge origin/master

1.Before working on the folder
git fetch
git merge

2.Do your work
git add .
Commit your changes

git commit -a
or
->git commit -a -m "This is the summary of what you did"

3.Push your changes to the repository in the server
git push origin master





#Pull changes from remote and merge
git pull 


Good source to check how to use git
https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud


