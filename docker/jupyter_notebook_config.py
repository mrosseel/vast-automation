## Hashed password to use for web authentication.
#
#  To generate, type in a python/IPython shell:
#
#    from notebook.auth import passwd; passwd()
#
#  The string should be of the form type:salt:hashed-password.
c.NotebookApp.password = u'sha1:8a9f684b1c03:9830877069eb6c4a32c24239f14c5cd58a51f7ba'
