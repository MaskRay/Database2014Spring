CENTER = new google.maps.LatLng(1.4634257573177956, 103.7522048939718)

class Map
  constructor = ->
    @ms = []
    @map = new google.maps.Map $('#map')[0],
      zoom: 12
      center: CENTER
      mapTypeId: google.maps.MapTypeId.ROADMAP
    @marker = new google.maps.Marker
      position: latlng
      icon: 'http://maps.google.com/mapfiles/marker.png'
    @marker.setMap map
    @map.setCenter marker.getPosition()

    google.maps.event.addListener map, 'click', (ev)->
      marker.position = ev.latLng
      marker.setMap map
      keypress()
      q = $('#query').val()
      lat = marker.position.lat()
      lng = marker.position.lng()
      $.getJSON '/search', {q, lat, lng}, (data)->
        console.log data
        for m in @ms
          m.setMap null
        @ms = []
        for p in data
          pos = new google.maps.LatLng p.lat, p.lng
          marker = new google.maps.Marker position: pos, icon: 'http://maps.google.com/mapfiles/marker_green.png'
          marker.setMap map
          marker.setTitle p.name
          @ms.push marker
