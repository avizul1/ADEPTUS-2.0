from flask.ext.wtf import Form
from wtforms import StringField, BooleanField, TextAreaField
from wtforms.validators import DataRequired

class LoginForm(Form):
    #if there is form with validadtor - to remove from comment state
    openid = TextAreaField('openid')#, validators=[DataRequired()]) 
    #remember_me = BooleanField('remember_me', default=False)
