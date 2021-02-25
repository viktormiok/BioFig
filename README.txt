#Clone the repository to your local computer as follows (substitute username for your username)
git clone username@idocdr1:/storage/cdr/Shared/BioFig.git BioFig.git

#check URLs of remote
git remote -v

#Create new files or modify them, then:
git add .

#Commit your changes
git commit -a
or
git commit -a -m "This is the summary of what you did"

#Push your changes to the repository in the server
git push origin master

#Pull changes safe way not merging
git fetch

#Pull changes from remote and merge
git pull 


Good source to check how to use git
https://www.atlassian.com/git/tutorials/learn-git-with-bitbucket-cloud