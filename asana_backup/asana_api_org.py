import asana
from asana.rest import ApiException
from pprint import pprint
import json
import urllib3
#urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

user_gid = "me"
workspace_gid = ''
team_gid = ''
cert_loc = r''

configuration = asana.Configuration()
#configuration.ssl_ca_cert = cert_loc
configuration.verify_ssl = False
configuration.access_token = ''
api_client = asana.ApiClient(configuration)


def get_user_info(api_client):
# create an instance of the API class
    users_api_instance = asana.UsersApi(api_client)
    user_gid = "me" # str | A string identifying a user. This can either be the string \"me\", an email, or the gid of a user.
    opts = {
        'opt_fields': "email,name,workspaces,workspaces.name,team", # list[str] | This endpoint returns a resource which excludes some properties by default. To include those optional properties, set this query parameter to a comma-separated list of the properties you wish to include.
    }

    try:
        # Get a user
        api_response = users_api_instance.get_user(user_gid, opts)
        pprint(api_response)
    except ApiException as e:
        print("Exception when calling UsersApi->get_user: %s\n" % e)
        
def get_workspace_info(api_client):
    workspaces_api_instance = asana.WorkspacesApi(api_client)
    workspace_gid = "" # str | Globally unique identifier for the workspace or organization.
    opts = {
        'opt_fields': "is_organization,name", 
    }

    try:
        # Get a workspace
        api_response = workspaces_api_instance.get_workspace(workspace_gid, opts)
        print("Workspace Info:")
        pprint(api_response)
    except ApiException as e:
        print("Exception when calling WorkspacesApi->get_workspace: %s\n" % e)
        
        
def get_user_team_info(api_client):
# create an instance of the API class
    team_memberships_api_instance = asana.TeamMembershipsApi(api_client)
    user_gid = "me" # str | A string identifying a user. This can either be the string \"me\", an email, or the gid of a user.
    workspace = workspace_gid # str | Globally unique identifier for the workspace.
    opts = {
        'limit': 50, # int | Results per page. The number of objects to return per page. The value must be between 1 and 100.
        'opt_fields': "is_admin,is_guest,is_limited_access,offset,path,team,team.name,uri,user,user.name", # list[str] | This endpoint returns a resource which excludes some properties by default. To include those optional properties, set this query parameter to a comma-separated list of the properties you wish to include.
    }

    try:
        # Get memberships from a user
        api_response = team_memberships_api_instance.get_team_memberships_for_user(user_gid, workspace, opts)
        for data in api_response:
            pprint(data)
    except ApiException as e:
        print("Exception when calling TeamMembershipsApi->get_team_memberships_for_user: %s\n" % e)

def get_single_project(api_client):
    
# create an instance of the API class
    projects_api_instance = asana.ProjectsApi(api_client)
    project_gid =  # str | Globally unique identifier for the project.
    opts = {
        'opt_fields': "archived,color,completed,completed_at,completed_by,completed_by.name,created_at,created_from_template,created_from_template.name,current_status,current_status.author,current_status.author.name,current_status.color,current_status.created_at,current_status.created_by,current_status.created_by.name,current_status.html_text,current_status.modified_at,current_status.text,current_status.title,current_status_update,current_status_update.resource_subtype,current_status_update.title,custom_field_settings,custom_field_settings.custom_field,custom_field_settings.custom_field.asana_created_field,custom_field_settings.custom_field.created_by,custom_field_settings.custom_field.created_by.name,custom_field_settings.custom_field.currency_code,custom_field_settings.custom_field.custom_label,custom_field_settings.custom_field.custom_label_position,custom_field_settings.custom_field.date_value,custom_field_settings.custom_field.date_value.date,custom_field_settings.custom_field.date_value.date_time,custom_field_settings.custom_field.default_access_level,custom_field_settings.custom_field.description,custom_field_settings.custom_field.display_value,custom_field_settings.custom_field.enabled,custom_field_settings.custom_field.enum_options,custom_field_settings.custom_field.enum_options.color,custom_field_settings.custom_field.enum_options.enabled,custom_field_settings.custom_field.enum_options.name,custom_field_settings.custom_field.enum_value,custom_field_settings.custom_field.enum_value.color,custom_field_settings.custom_field.enum_value.enabled,custom_field_settings.custom_field.enum_value.name,custom_field_settings.custom_field.format,custom_field_settings.custom_field.has_notifications_enabled,custom_field_settings.custom_field.id_prefix,custom_field_settings.custom_field.is_formula_field,custom_field_settings.custom_field.is_global_to_workspace,custom_field_settings.custom_field.is_value_read_only,custom_field_settings.custom_field.multi_enum_values,custom_field_settings.custom_field.multi_enum_values.color,custom_field_settings.custom_field.multi_enum_values.enabled,custom_field_settings.custom_field.multi_enum_values.name,custom_field_settings.custom_field.name,custom_field_settings.custom_field.number_value,custom_field_settings.custom_field.people_value,custom_field_settings.custom_field.people_value.name,custom_field_settings.custom_field.precision,custom_field_settings.custom_field.privacy_setting,custom_field_settings.custom_field.representation_type,custom_field_settings.custom_field.resource_subtype,custom_field_settings.custom_field.text_value,custom_field_settings.custom_field.type,custom_field_settings.is_important,custom_field_settings.parent,custom_field_settings.parent.name,custom_field_settings.project,custom_field_settings.project.name,custom_fields,custom_fields.date_value,custom_fields.date_value.date,custom_fields.date_value.date_time,custom_fields.display_value,custom_fields.enabled,custom_fields.enum_options,custom_fields.enum_options.color,custom_fields.enum_options.enabled,custom_fields.enum_options.name,custom_fields.enum_value,custom_fields.enum_value.color,custom_fields.enum_value.enabled,custom_fields.enum_value.name,custom_fields.id_prefix,custom_fields.is_formula_field,custom_fields.multi_enum_values,custom_fields.multi_enum_values.color,custom_fields.multi_enum_values.enabled,custom_fields.multi_enum_values.name,custom_fields.name,custom_fields.number_value,custom_fields.representation_type,custom_fields.text_value,custom_fields.type,default_access_level,default_view,due_date,due_on,followers,followers.name,html_notes,icon,members,members.name,minimum_access_level_for_customization,minimum_access_level_for_sharing,modified_at,name,notes,owner,permalink_url,privacy_setting,project_brief,public,start_on,team,team.name,workspace,workspace.name", # list[str] | This endpoint returns a resource which excludes some properties by default. To include those optional properties, set this query parameter to a comma-separated list of the properties you wish to include.
}

    try:
        # Get a project
        api_response = projects_api_instance.get_project(project_gid, opts)
        
        with open("projects.txt", 'w') as f:
            pprint(api_response, stream=f)
        
        #pprint(api_response)
        
        
        
    except ApiException as e:
        print("Exception when calling ProjectsApi->get_project: %s\n" % e)
        
    

def get_projects(api_client):
    projects_api_instance = asana.ProjectsApi(api_client)
    team_gid = "" # str | Globally unique identifier for the team.
    opts = {
        'limit': 100, # int | Results per page. The number of objects to return per page. The value must be between 1 and 100.
       
        'archived': True, # bool | Only return projects whose `archived` field takes on the value of this parameter.
        'opt_fields': "archived,color,completed,completed_at,completed_by,completed_by.name,created_at,created_from_template,created_from_template.name,current_status,current_status.author,current_status.author.name,current_status.color,current_status.created_at,current_status.created_by,current_status.created_by.name,current_status.html_text,current_status.modified_at,current_status.text,current_status.title,current_status_update,current_status_update.resource_subtype,current_status_update.title,custom_field_settings,custom_field_settings.custom_field,custom_field_settings.custom_field.asana_created_field,custom_field_settings.custom_field.created_by,custom_field_settings.custom_field.created_by.name,custom_field_settings.custom_field.currency_code,custom_field_settings.custom_field.custom_label,custom_field_settings.custom_field.custom_label_position,custom_field_settings.custom_field.date_value,custom_field_settings.custom_field.date_value.date,custom_field_settings.custom_field.date_value.date_time,custom_field_settings.custom_field.default_access_level,custom_field_settings.custom_field.description,custom_field_settings.custom_field.display_value,custom_field_settings.custom_field.enabled,custom_field_settings.custom_field.enum_options,custom_field_settings.custom_field.enum_options.color,custom_field_settings.custom_field.enum_options.enabled,custom_field_settings.custom_field.enum_options.name,custom_field_settings.custom_field.enum_value,custom_field_settings.custom_field.enum_value.color,custom_field_settings.custom_field.enum_value.enabled,custom_field_settings.custom_field.enum_value.name,custom_field_settings.custom_field.format,custom_field_settings.custom_field.has_notifications_enabled,custom_field_settings.custom_field.id_prefix,custom_field_settings.custom_field.is_formula_field,custom_field_settings.custom_field.is_global_to_workspace,custom_field_settings.custom_field.is_value_read_only,custom_field_settings.custom_field.multi_enum_values,custom_field_settings.custom_field.multi_enum_values.color,custom_field_settings.custom_field.multi_enum_values.enabled,custom_field_settings.custom_field.multi_enum_values.name,custom_field_settings.custom_field.name,custom_field_settings.custom_field.number_value,custom_field_settings.custom_field.people_value,custom_field_settings.custom_field.people_value.name,custom_field_settings.custom_field.precision,custom_field_settings.custom_field.privacy_setting,custom_field_settings.custom_field.representation_type,custom_field_settings.custom_field.resource_subtype,custom_field_settings.custom_field.text_value,custom_field_settings.custom_field.type,custom_field_settings.is_important,custom_field_settings.parent,custom_field_settings.parent.name,custom_field_settings.project,custom_field_settings.project.name,custom_fields,custom_fields.date_value,custom_fields.date_value.date,custom_fields.date_value.date_time,custom_fields.display_value,custom_fields.enabled,custom_fields.enum_options,custom_fields.enum_options.color,custom_fields.enum_options.enabled,custom_fields.enum_options.name,custom_fields.enum_value,custom_fields.enum_value.color,custom_fields.enum_value.enabled,custom_fields.enum_value.name,custom_fields.id_prefix,custom_fields.is_formula_field,custom_fields.multi_enum_values,custom_fields.multi_enum_values.color,custom_fields.multi_enum_values.enabled,custom_fields.multi_enum_values.name,custom_fields.name,custom_fields.number_value,custom_fields.representation_type,custom_fields.text_value,custom_fields.type,default_access_level,default_view,due_date,due_on,followers,followers.name,html_notes,icon,members,members.name,minimum_access_level_for_customization,minimum_access_level_for_sharing,modified_at,name,notes,offset,owner,path,permalink_url,privacy_setting,project_brief,public,start_on,team,team.name,uri,workspace,workspace.name", # list[str] | This endpoint returns a resource which excludes some properties by default. To include those optional properties, set this query parameter to a comma-separated list of the properties you wish to include.
}

    try:
        # Get a team's projects
        api_response = projects_api_instance.get_projects_for_team(team_gid, opts)
        
        
        project_list = list(api_response)
        
        with open('project_list_archived.json','w') as f:
            json.dump(project_list, f, indent=4)
        
        '''
        with open("projects.txt", 'w') as f:
            for data in list(api_response):
                pprint(data, stream=f)
        '''
    except ApiException as e:
        print("Exception when calling ProjectsApi->get_projects_for_team: %s\n" % e)
        
    return api_response
    
#get_user_info(api_client)
#get_workspace_info(api_client)
#get_teams_in_workspace_info(api_client)
#projects = get_user_team_info(api_client)
#project = get_single_project(api_client)
projects = get_projects(api_client)
