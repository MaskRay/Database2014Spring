#!/usr/bin/env ruby

require 'open-uri'
require 'socket'
require 'cgi'
require 'json'
begin
  require 'sass'
  require 'slim'
  require 'coffee-script'
  require 'sinatra'
  require 'sinatra/reloader'
rescue LoadError
  STDERR.puts 'gem install sinatra sinatra-contrib sass slim coffee-script'
end

# Main

configure :development do
  register Sinatra::Reloader
end

set :static, true
set :public_folder, File.dirname(__FILE__)
set :views, File.dirname(__FILE__)
set :bind, '0'
set :port, 4567

## Routes

get '/css/*.css' do
  sass params[:splat][0].to_sym
end

get '/js/*.js' do
  coffee params[:splat][0].to_sym
end

get '/' do
  slim :index
end

get '/main' do
  slim :main
end

get '/report' do
  slim :report
end

get '/autocomplete' do
  content_type :json

  query = CGI::parse request.query_string
  q = query['query'][0].strip
  s = TCPSocket.new '0', 1337
  s.write "1\n#{q}\n"
  s.close_write
  r = s.read
  print r
  {suggestions: JSON::parse(r).map {|x| x['name'] }}.to_json
end

get '/search' do
  content_type :json
  query = CGI::parse request.query_string
  q = query.key?('q') ? query['q'][0].strip : ''
  puts "q: #{q}"
  s = TCPSocket.new '0', 1337
  if query.key? 'lat'
    lat = query['lat'][0]
    lng = query['lng'][0]
    s.write "0\n#{lat} #{lng} #{q}\n"
  else
    s.write "1\n#{q}\n"
  end
  s.close_write

  r = s.read
  print r
  r
end
