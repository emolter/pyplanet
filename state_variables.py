from __future__ import print_function, absolute_import, division


def init_state_variables(mode, **kwargs):
    state_vars = {'batch_mode': False,
                  'write_output_files': True,
                  'write_log_file': True,
                  'plot': True,
                  'verbose': False,
                  'super_quiet': False,
                  'generate_alpha': False,
                  'use_existing_alpha': False,
                  'scale_existing_alpha': False,
                  'scale_file_name': None,
                  'output_type': 'frequency'
                  }

    if mode == 'batch':
        state_vars['batch_mode'] = True
        state_vars['self.plot'] = False
        state_vars['verbose'] = False
        state_vars['super_quiet'] = True
        state_vars['write_log_file'] = False
    elif mode == 'mcmc':
        state_vars['plot'] = False
        state_vars['verbose'] = False
        state_vars['super_quiet'] = True
        state_vars['write_log_file'] = False
        state_vars['write_output_files'] = False
        state_vars['scale_existing_alpha'] = True
        state_vars['scale_file_name'] = 'Scratch/scale.dat'
    elif mode == 'use_alpha':
        state_vars['plot'] = False
        state_vars['verbose'] = False
        state_vars['use_existing_alpha'] = True
    elif mode == 'scale_alpha':
        state_vars['plot'] = False
        state_vars['verbose'] = False
        state_vars['scale_existing_alpha'] = True
        state_vars['scale_file_name'] = 'Scratch/scale.dat'

    # update based on provided kwargs
    for k in kwargs:
        if k in state_vars.keys():
            state_vars[k] = kwargs[k]
        else:
            print("{} keyword not found.".format(k))
            raise ValueError("Aborting since you probably wanted this keyword")

    if state_vars['super_quiet']:
        state_vars['verbose'] = False

    return state_vars
