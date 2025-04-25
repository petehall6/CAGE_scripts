import asana
from asana.rest import ApiException
import pprint
# Configure personal access token
configuration = asana.Configuration()
configuration.verify_ssl = False
configuration.access_token = '2/1203461498038908/1209949224473555:b678631bb89636dc50465990550229e4'
api_client = asana.ApiClient(configuration)


#get user info

def get_user_info(api_client):
    
    users_api_instance = asana.UsersApi(api_client)

    user_gid = "me" # str | A string identifying a user. This can either be the string \"me\", an email, or the gid of a user.
    opts = {
        'opt_fields': "name,workspaces,workspaces.name",
    }

    try:
        # Get a user
        api_response = users_api_instance.get_user(user_gid, opts)
        pprint(api_response)
    except ApiException as e:
        print("Exception when calling UsersApi->get_user: %s\n" % e)


get_user_info(api_client)