module.exports = (grunt) ->
  'use strict'

  grunt.initConfig
    pkg: grunt.file.readJSON('package.json')

    copy:
      dev:
        files: [ {
          expand: true
          src: ['blockly/block*_compressed.js', 'blockly/msg/js/en.js', 'blockly/blocks/*', 'assets/*']
          dest: 'build'
        } ]

    coffee:
      dev:
        files: [
          expand: true
          cwd: 'js'
          src: ['**/*.coffee', '!**/_*.coffee']
          dest: 'build/js'
          ext: '.js'
        ]

    sass:
      dev:
        options:
          loadPath: 'bower_components/foundation/scss'
        files: [
          expand: true
          cwd: 'css'
          src: ['**/*.sass', '!**/_*.sass']
          dest: 'build/css'
          ext: '.css'
        ]

    slim:
      dev:
        files: [ {
          expand: true
          src: ['*.slim', '!_*.slim']
          dest: 'build'
          ext: '.html'
        } ]

    watch:
      options:
        livereload: true
        spawn: false

      coffee:
        files: ['js/**/*.coffee']
        tasks: ['coffee:dev']

      sass:
        files: ['css/**/*.sass']
        tasks: ['sass:dev']

      slim:
        files: ['*.slim']
        tasks: ['slim:dev']

    connect:
      dev:
        options:
          port: 9999
          hostname: '*'
          base: 'build'
          livereload: true
          open: true

    clean:
      dev: ['build/*']

    bower:
      dev:
        options:
          targetDir: '.'
          layout: 'byType'

  grunt.loadNpmTasks 'grunt-bower-task'
  grunt.loadNpmTasks 'grunt-contrib-coffee'
  grunt.loadNpmTasks 'grunt-contrib-connect'
  grunt.loadNpmTasks 'grunt-contrib-copy'
  grunt.loadNpmTasks 'grunt-contrib-sass'
  grunt.loadNpmTasks 'grunt-contrib-watch'
  grunt.loadNpmTasks 'grunt-shell'
  grunt.loadNpmTasks 'grunt-slim'

  grunt.registerTask 'dev', ['copy:dev', 'coffee:dev', 'sass:dev', 'slim:dev']
  grunt.registerTask 'default', ['dev', 'connect', 'watch']
