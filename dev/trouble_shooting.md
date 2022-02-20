## ドキュメント公開の流れ
- docsに移動
- `make html`を実行
- Readthedocsに移動し，ビルド

## pypi公開
- 

## Q & A
```
Traceback (most recent call last):
  File "/Users/ogyogugige/anaconda3/bin/sphinx-build", line 5, in <module>
    from sphinx.cmd.build import main
  File "/Users/ogyogugige/anaconda3/lib/python3.6/site-packages/sphinx/cmd/build.py", line 25, in <module>
    from sphinx.application import Sphinx
  File "/Users/ogyogugige/anaconda3/lib/python3.6/site-packages/sphinx/application.py", line 34, in <module>
    from sphinx.domains import Domain, Index
  File "/Users/ogyogugige/anaconda3/lib/python3.6/site-packages/sphinx/domains/__init__.py", line 24, in <module>
    from sphinx.roles import XRefRole
  File "/Users/ogyogugige/anaconda3/lib/python3.6/site-packages/sphinx/roles.py", line 20, in <module>
    from sphinx.util.docutils import ReferenceRole, SphinxRole
  File "/Users/ogyogugige/anaconda3/lib/python3.6/site-packages/sphinx/util/docutils.py", line 44, in <module>
    __version_info__ = version.parse(docutils.__version__).release
AttributeError: 'Version' object has no attribute 'release'
```
以下を実行することで，解消可能
```
pip install -U packaging
```
