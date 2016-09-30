## Simple instructions to build LuaJIT 2.1
##
## NOTES: Post install you may want to delete the dynamic libraries,
## so that the static library is linked to gkyl.
##
## If "git checkout v2.1' does not work just download the ZIP file for
## LuaJIT 2.1beta1 and install that instead.

git clone http://luajit.org/git/luajit-2.0.git
cd luajit-2.0
git checkout v2.1
make
make install PREFIX=$HOME/gkylsoft/luajit