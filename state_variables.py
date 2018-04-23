from __future__ import print_function, absolute_import, division


def init_state_variables(mode, **kwargs):
    state_vars = {'batch_mode': False,
                  'write_output_files': True,
                  'write_log_file': True,
                  'plot': True,
                  'verbose': False,
                  'generate_alpha': False,
                  'use_alpha': False,
                  'output_type': 'frequency'
                  }

    if mode == 'batch':
        state_vars['batch_mode'] = True
        state_vars['self.plot'] = False
        state_vars['verbose'] = False
        state_vars['write_log_file'] = False
    elif mode == 'mcmc':
        state_vars['plot'] = False
        state_vars['verbose'] = False
        state_vars['write_log_file'] = False
        state_vars['write_output_files'] = False
        state_vars['use_alpha'] = True
    elif mode == 'generate_alpha':
        state_vars['write_output_files'] = False
        state_vars['generate_alpha'] = True

    # update based on provided kwargs
    for k in kwargs:
        if k in state_vars.keys():
            state_vars[k] = kwargs[k]

    return state_vars
